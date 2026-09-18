#!/usr/bin/env python3
"""Extract a gene by its GFF3 ID and write all transcripts as BED12."""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


TRANSCRIPT_TYPES = {
    "mrna",
    "transcript",
    "ncrna",
    "lncrna",
    "mirna",
    "rrna",
    "snrna",
    "snorna",
    "trna",
}
@dataclass(frozen=True)
class Feature:
    seqid: str
    kind: str
    start: int  # GFF3: 1-based, inclusive
    end: int
    strand: str
    attributes: dict[str, str]

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> list[str]:
        value = self.attributes.get("Parent", "")
        return [item for item in value.split(",") if item]


def open_text(path: Path):
    """Open plain-text or gzip-compressed GFF3."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open("r", encoding="utf-8-sig")


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in text.strip().strip(";").split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        else:
            parts = field.strip().split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
            value = value.strip('"')
        attributes[unquote(key.strip())] = unquote(value.strip())
    return attributes


def read_gff3(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"第 {line_number} 行不是合法的 9 列 GFF3")
            try:
                start, end = int(columns[3]), int(columns[4])
            except ValueError as exc:
                raise ValueError(f"第 {line_number} 行的坐标不是整数") from exc
            if start < 1 or end < start:
                raise ValueError(f"第 {line_number} 行的坐标范围无效")
            features.append(
                Feature(
                    seqid=columns[0],
                    kind=columns[2],
                    start=start,
                    end=end,
                    strand=columns[6] if columns[6] in {"+", "-"} else ".",
                    attributes=parse_attributes(columns[8]),
                )
            )
    return features


def matches_gene_id(feature: Feature, gene_id: str) -> bool:
    """Match only the exact ID of a GFF3 gene feature."""
    return feature.kind.lower() == "gene" and feature.feature_id == gene_id


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def safe_bed_name(value: str) -> str:
    return value.replace("\t", "_").replace(" ", "_")


def make_bed12(features: list[Feature], gene_id_query: str) -> list[str]:
    genes = [feature for feature in features if matches_gene_id(feature, gene_id_query)]
    if not genes:
        raise LookupError(
            f"没有找到 ID={gene_id_query} 的 gene；请使用 GFF3 gene 行中 ID= 后的完整值"
        )
    if len(genes) > 1:
        raise LookupError(f"GFF3 中存在 {len(genes)} 个 ID={gene_id_query} 的 gene，ID 不唯一")

    gene = genes[0]
    gene_id = gene.feature_id
    if not gene_id:
        raise ValueError("匹配到的 gene feature 没有 ID 属性，无法追踪其子级")

    transcripts = {
        feature.feature_id: feature
        for feature in features
        if feature.kind.lower() in TRANSCRIPT_TYPES
        and gene_id in feature.parents
        and feature.feature_id
    }

    # Some annotations attach exon/CDS directly to the gene rather than a transcript.
    if not transcripts:
        transcripts = {gene_id: gene}

    children: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        for parent in feature.parents:
            if parent in transcripts:
                children[parent].append(feature)

    gene_name = next(
        (gene.attributes[key] for key in ("Name", "gene_name", "locus_tag") if key in gene.attributes),
        gene_id,
    )
    lines: list[str] = []
    for transcript_id, transcript in transcripts.items():
        transcript_children = children.get(transcript_id, [])
        exons = [f for f in transcript_children if f.kind.lower() == "exon"]
        cds = [f for f in transcript_children if f.kind.lower() == "cds"]
        if not exons:
            exons = cds
        if not exons:
            print(f"警告：跳过没有 exon/CDS 的转录本 {transcript_id}", file=sys.stderr)
            continue
        if any(f.seqid != transcript.seqid for f in exons):
            raise ValueError(f"转录本 {transcript_id} 的外显子分布在不同序列上")

        # Convert GFF3 1-based inclusive coordinates to BED 0-based half-open.
        blocks = merge_intervals([(f.start - 1, f.end) for f in exons])
        chrom_start, chrom_end = blocks[0][0], blocks[-1][1]
        cds_on_chrom = [f for f in cds if f.seqid == transcript.seqid]
        if cds_on_chrom:
            thick_start = min(f.start for f in cds_on_chrom) - 1
            thick_end = max(f.end for f in cds_on_chrom)
        else:
            thick_start = thick_end = chrom_start

        transcript_name = next(
            (
                transcript.attributes[key]
                for key in ("Name", "transcript_name")
                if key in transcript.attributes
            ),
            transcript_id,
        )
        name = safe_bed_name(f"{gene_name}|{transcript_name}")
        block_sizes = ",".join(str(end - start) for start, end in blocks) + ","
        block_starts = ",".join(str(start - chrom_start) for start, _ in blocks) + ","
        lines.append(
            "\t".join(
                map(
                    str,
                    (
                        transcript.seqid,
                        chrom_start,
                        chrom_end,
                        name,
                        0,
                        transcript.strand,
                        thick_start,
                        thick_end,
                        "0,102,204",
                        len(blocks),
                        block_sizes,
                        block_starts,
                    ),
                )
            )
        )

    if not lines:
        raise LookupError(f"基因 {gene_id_query} 没有可输出的 exon 或 CDS")
    return lines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="按 gene ID 从 GFF3 提取全部转录本，生成 pyGenomeTracks 用 BED12。"
    )
    parser.add_argument("gff3", type=Path, help="输入的 .gff3 或 .gff3.gz 文件")
    parser.add_argument("gene_id", help="GFF3 gene 行中 ID= 后的完整值（精确匹配）")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("genes.bed"), help="输出文件（默认：genes.bed）"
    )
    parser.add_argument("--force", action="store_true", help="覆盖已经存在的输出文件")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.gff3.is_file():
        print(f"错误：找不到 GFF3 文件：{args.gff3}", file=sys.stderr)
        return 1
    if args.gff3.resolve() == args.output.resolve():
        print("错误：输出文件不能与输入 GFF3 是同一个文件", file=sys.stderr)
        return 1
    if args.output.exists() and not args.force:
        print(f"错误：输出文件已存在：{args.output}（使用 --force 覆盖）", file=sys.stderr)
        return 1
    try:
        lines = make_bed12(read_gff3(args.gff3), args.gene_id)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    except (OSError, ValueError, LookupError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    print(f"已输出 {len(lines)} 个转录本：{args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
