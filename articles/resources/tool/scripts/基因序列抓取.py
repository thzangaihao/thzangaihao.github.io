#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from urllib.parse import unquote

"""
特定基因组序列精确提取与蛋白翻译工具

功能：
1. 支持手动输入，或从 csv/tsv/txt 文件批量读取基因 ID。
2. 基于 GFF/GTF 注释提取基因全长 DNA 序列。
3. 基于每个转录本的 CDS 拼接并翻译蛋白序列。
4. 蛋白输出支持：所有转录本、最长转录本、两者都输出。
5. FASTA 输出为两行格式，方便后续脚本处理。
"""


TRANSCRIPT_FEATURES = {"mrna", "transcript", "rna"}


def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))


def check_env():
    if not shutil.which("samtools"):
        print("错误: 未找到 samtools。请先安装，例如: conda install -c bioconda samtools")
        sys.exit(1)


def interactive_select(files, desc, display_root):
    if not files:
        print(f"未找到任何 {desc}。")
        return []

    files = sorted(files)
    print(f"\n在当前目录及子目录中扫描到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        rel_path = os.path.relpath(f, display_root)
        print(f"  [{i}] {rel_path}")

    while True:
        choice = input(
            f"\n请选择 {desc} (单选:1, 多选:1,3, 范围:1-3, 全部:all, 退出:q): "
        ).strip().lower()
        if choice == "q":
            sys.exit(0)
        if choice == "all":
            return files

        try:
            selected = []
            parts = choice.replace(" ", "").split(",")
            for part in parts:
                if "-" in part:
                    start, end = map(int, part.split("-", 1))
                    selected.extend(files[start - 1:end])
                else:
                    selected.append(files[int(part) - 1])
            return list(dict.fromkeys(selected))
        except (ValueError, IndexError):
            print("输入格式错误，请重新选择。")


def clean_feature_id(raw_id):
    if not raw_id:
        return ""
    cleaned = unquote(str(raw_id)).strip().strip('"').strip("'")
    cleaned = cleaned.split("|")[-1]
    cleaned = re.sub(
        r"^(gene-|rna-|mrna\.|mrna-|transcript-|cds-|id:)",
        "",
        cleaned,
        flags=re.IGNORECASE,
    )
    return cleaned


def split_attr_values(value):
    if not value:
        return []
    values = []
    for item in re.split(r"[,;]", value):
        item = clean_feature_id(item)
        if item:
            values.append(item)
    return values


def parse_attributes(attr_text):
    attrs = {}
    for item in attr_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"').strip("'")
    return attrs


def collect_aliases(attrs, keys):
    aliases = set()
    for key in keys:
        aliases.update(split_attr_values(attrs.get(key)))
    return aliases


def reverse_complement(seq):
    trans = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return seq.translate(trans)[::-1]


def translate_dna_to_protein(seq):
    codon_table = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
    }
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i + 3]
        protein.append(codon_table.get(codon, "X"))
    return "".join(protein)


def find_gene_input_files(base_dir):
    patterns = ["*.txt", "*.csv", "*.tsv"]
    files = []
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(base_dir, "**", pattern), recursive=True))
    return sorted(set(files))


def read_gene_ids_from_file(path):
    ext = os.path.splitext(path)[1].lower()
    genes = []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        if ext == ".csv":
            reader = csv.reader(handle)
            for row in reader:
                genes.extend(row)
        elif ext == ".tsv":
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                genes.extend(row)
        else:
            text = handle.read()
            genes.extend(re.split(r"[\s,;，；]+", text))

    cleaned = [clean_feature_id(g) for g in genes if clean_feature_id(g)]
    return list(dict.fromkeys(cleaned))


def get_target_genes(base_dir):
    print("\n--- 输入目标基因 ---")
    print("  [1] 手动输入基因 ID，多个 ID 用逗号分隔")
    print("  [2] 从 csv/tsv/txt 文件批量读取")
    print("  [3] 提取注释文件中的所有基因")
    input_mode = input("请选择输入方式 (1/2/3) [默认 3]: ").strip()
    if input_mode not in {"1", "2", "3"}:
        input_mode = "3"

    if input_mode == "3":
        targets = None
        print("已选择: 提取注释文件中的所有基因。")
    elif input_mode == "1":
        genes_input = input("请输入基因 ID: ").strip()
        targets = [clean_feature_id(g) for g in re.split(r"[,，\s;；]+", genes_input) if g.strip()]
    else:
        input_files = find_gene_input_files(base_dir)
        selected_files = interactive_select(input_files, "基因列表文件 (.csv/.tsv/.txt)", base_dir)
        targets = []
        for path in selected_files:
            targets.extend(read_gene_ids_from_file(path))
        targets = list(dict.fromkeys(targets))

    if targets is not None and not targets:
        print("未读取到任何基因 ID，程序退出。")
        sys.exit(0)

    if targets is not None:
        print(f"已读取 {len(targets)} 个目标 ID。")

    print("\n--- 选择输出类型 ---")
    print("  [1] 仅输出基因全长核苷酸序列 (DNA)")
    print("  [2] 仅输出 CDS 翻译的蛋白序列 (Protein)")
    print("  [3] 两者都要")
    output_mode = input("请选择 (1/2/3) [默认 2]: ").strip()
    if output_mode not in {"1", "2", "3"}:
        output_mode = "2"

    protein_mode = None
    if output_mode in {"2", "3"}:
        print("\n--- 选择蛋白转录本范围 ---")
        print("  [1] 所有转录本")
        print("  [2] 每个基因的最长转录本")
        print("  [3] 两者都要")
        protein_mode = input("请选择 (1/2/3) [默认 3]: ").strip()
        if protein_mode not in {"1", "2", "3"}:
            protein_mode = "3"

    return (None if targets is None else set(targets)), output_mode, protein_mode


def parse_annotation(gff_file):
    genes = []
    transcripts = []
    cds_records = []

    with open(gff_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attr_text = parts
            feature_lower = feature.lower()
            attrs = parse_attributes(attr_text)
            start_i, end_i = int(start), int(end)

            if feature_lower == "gene":
                aliases = collect_aliases(
                    attrs,
                    ["ID", "Name", "gene", "gene_id", "gene_name", "locus_tag", "Alias"],
                )
                gene_id = clean_feature_id(attrs.get("ID")) or next(iter(aliases), "")
                if gene_id:
                    aliases.add(gene_id)
                genes.append({
                    "id": gene_id,
                    "aliases": aliases,
                    "chrom": chrom,
                    "strand": strand,
                    "coord": (start_i, end_i),
                })
            elif feature_lower in TRANSCRIPT_FEATURES:
                aliases = collect_aliases(
                    attrs,
                    ["ID", "Name", "transcript_id", "transcript_name", "Alias"],
                )
                parent_aliases = collect_aliases(
                    attrs,
                    ["Parent", "gene", "gene_id", "gene_name", "locus_tag"],
                )
                transcript_id = clean_feature_id(attrs.get("ID")) or clean_feature_id(attrs.get("transcript_id"))
                if transcript_id:
                    aliases.add(transcript_id)
                transcripts.append({
                    "id": transcript_id,
                    "aliases": aliases,
                    "parent_aliases": parent_aliases,
                    "chrom": chrom,
                    "strand": strand,
                    "coord": (start_i, end_i),
                })
            elif feature_lower == "cds":
                parent_aliases = collect_aliases(
                    attrs,
                    ["Parent", "transcript_id", "gene", "gene_id", "gene_name", "locus_tag"],
                )
                cds_records.append({
                    "parents": parent_aliases,
                    "chrom": chrom,
                    "strand": strand,
                    "coord": (start_i, end_i),
                    "phase": phase,
                })

    return genes, transcripts, cds_records


def choose_record_id(primary_id, aliases, fallback):
    if primary_id:
        return primary_id
    aliases = sorted(alias for alias in aliases if alias)
    return aliases[0] if aliases else fallback


def finalize_found_genes(found):
    valid = {}
    for target, info in found.items():
        cds_coords = [coord for tx in info["transcripts"].values() for coord in tx["cds_coords"]]
        if not info.get("chrom") and info["transcripts"]:
            first_tx = next(iter(info["transcripts"].values()))
            info["chrom"] = first_tx["chrom"]
            info["strand"] = first_tx["strand"]
        if not info.get("gene_coord") and cds_coords:
            info["gene_coord"] = (min(c[0] for c in cds_coords), max(c[1] for c in cds_coords))
        if info.get("chrom") or cds_coords:
            valid[target] = info
    return valid


def parse_gff_for_all_genes(gff_file):
    print("\n正在解析注释文件并定位所有基因/转录本...")
    genes, transcripts, cds_records = parse_annotation(gff_file)

    found = {}
    gene_alias_to_id = {}
    transcript_owner = {}

    for index, gene in enumerate(genes, 1):
        gene_id = choose_record_id(gene["id"], gene["aliases"], f"gene_{index}")
        found.setdefault(gene_id, {
            "chrom": gene["chrom"],
            "strand": gene["strand"],
            "gene_coord": gene["coord"],
            "transcripts": {},
        })
        aliases = set(gene["aliases"])
        aliases.add(gene_id)
        for alias in aliases:
            if alias:
                gene_alias_to_id.setdefault(alias, gene_id)

    for index, tx in enumerate(transcripts, 1):
        tx_aliases = set(tx["aliases"])
        if tx["id"]:
            tx_aliases.add(tx["id"])

        gene_id = None
        for parent in tx["parent_aliases"]:
            if parent in gene_alias_to_id:
                gene_id = gene_alias_to_id[parent]
                break
        if gene_id is None:
            gene_id = choose_record_id("", tx["parent_aliases"], f"gene_from_transcript_{index}")

        tx_id = choose_record_id(tx["id"], tx_aliases, f"{gene_id}.transcript_{index}")
        found.setdefault(gene_id, {
            "chrom": tx["chrom"],
            "strand": tx["strand"],
            "gene_coord": tx["coord"],
            "transcripts": {},
        })
        if not found[gene_id].get("gene_coord"):
            found[gene_id]["gene_coord"] = tx["coord"]
        found[gene_id]["transcripts"].setdefault(tx_id, {
            "chrom": tx["chrom"],
            "strand": tx["strand"],
            "coord": tx["coord"],
            "cds_coords": [],
        })

        for alias in tx_aliases:
            if alias:
                transcript_owner[alias] = (gene_id, tx_id)
        for parent in tx["parent_aliases"]:
            if parent:
                gene_alias_to_id.setdefault(parent, gene_id)

    for index, cds in enumerate(cds_records, 1):
        owners = []
        for parent in cds["parents"]:
            if parent in transcript_owner:
                owners.append(transcript_owner[parent])

        if not owners:
            for parent in cds["parents"]:
                if parent in gene_alias_to_id:
                    gene_id = gene_alias_to_id[parent]
                    owners.append((gene_id, parent))

        if not owners:
            gene_id = choose_record_id("", cds["parents"], f"cds_parent_{index}")
            owners.append((gene_id, gene_id))

        for gene_id, tx_id in dict.fromkeys(owners):
            found.setdefault(gene_id, {
                "chrom": cds["chrom"],
                "strand": cds["strand"],
                "gene_coord": None,
                "transcripts": {},
            })
            found[gene_id]["transcripts"].setdefault(tx_id, {
                "chrom": cds["chrom"],
                "strand": cds["strand"],
                "coord": None,
                "cds_coords": [],
            })
            found[gene_id]["transcripts"][tx_id]["cds_coords"].append(cds["coord"])

    valid = finalize_found_genes(found)
    print(f"在注释文件中定位到 {len(valid)} 个可提取的基因/基因样记录。")
    return valid


def parse_gff_for_targets(gff_file, target_genes):
    if target_genes is None:
        return parse_gff_for_all_genes(gff_file)

    print("\n正在解析注释文件并定位目标基因/转录本...")
    genes, transcripts, cds_records = parse_annotation(gff_file)

    target_to_gene_aliases = defaultdict(set)
    found = {}
    transcript_owner = {}

    for gene in genes:
        matched_targets = gene["aliases"] & target_genes
        for target in matched_targets:
            target_to_gene_aliases[target].update(gene["aliases"])
            found.setdefault(target, {
                "chrom": gene["chrom"],
                "strand": gene["strand"],
                "gene_coord": gene["coord"],
                "transcripts": {},
            })

    for tx in transcripts:
        tx_aliases = set(tx["aliases"])
        tx_aliases.add(tx["id"])
        matched_targets = tx_aliases & target_genes
        matched_targets.update(tx["parent_aliases"] & target_genes)

        for target, gene_aliases in target_to_gene_aliases.items():
            if tx["parent_aliases"] & gene_aliases:
                matched_targets.add(target)

        for target in matched_targets:
            tx_id = tx["id"] or next(iter(tx_aliases - {""}), target)
            found.setdefault(target, {
                "chrom": tx["chrom"],
                "strand": tx["strand"],
                "gene_coord": tx["coord"],
                "transcripts": {},
            })
            if not found[target].get("gene_coord"):
                found[target]["gene_coord"] = tx["coord"]
            found[target]["transcripts"].setdefault(tx_id, {
                "chrom": tx["chrom"],
                "strand": tx["strand"],
                "coord": tx["coord"],
                "cds_coords": [],
            })
            for alias in tx_aliases:
                if alias:
                    transcript_owner[alias] = (target, tx_id)

    for target in target_genes:
        transcript_owner.setdefault(target, (target, target))

    for cds in cds_records:
        attached = False
        for parent in cds["parents"]:
            if parent in transcript_owner:
                target, tx_id = transcript_owner[parent]
                found.setdefault(target, {
                    "chrom": cds["chrom"],
                    "strand": cds["strand"],
                    "gene_coord": None,
                    "transcripts": {},
                })
                found[target]["transcripts"].setdefault(tx_id, {
                    "chrom": cds["chrom"],
                    "strand": cds["strand"],
                    "coord": None,
                    "cds_coords": [],
                })
                found[target]["transcripts"][tx_id]["cds_coords"].append(cds["coord"])
                attached = True

        if attached:
            continue

        direct_targets = cds["parents"] & target_genes
        for target in direct_targets:
            found.setdefault(target, {
                "chrom": cds["chrom"],
                "strand": cds["strand"],
                "gene_coord": None,
                "transcripts": {},
            })
            found[target]["transcripts"].setdefault(target, {
                "chrom": cds["chrom"],
                "strand": cds["strand"],
                "coord": None,
                "cds_coords": [],
            })
            found[target]["transcripts"][target]["cds_coords"].append(cds["coord"])

    valid = finalize_found_genes(found)

    print(f"共输入 {len(target_genes)} 个目标 ID，在注释文件中定位到 {len(valid)} 个。")
    missing = target_genes - set(valid.keys())
    if missing:
        print("以下目标未在注释文件中找到: " + ", ".join(sorted(missing)))
    return valid


def extract_sequence_chunk(fasta_file, chrom, start, end):
    result = subprocess.run(
        ["samtools", "faidx", fasta_file, f"{chrom}:{start}-{end}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    lines = result.stdout.strip().splitlines()
    if len(lines) > 1:
        return "".join(lines[1:])
    return ""


def format_fasta(header, seq):
    return f"{header}\n{seq}\n"


def build_protein_records(fasta_file, gene_id, info):
    records = []
    for tx_id, tx in sorted(info["transcripts"].items()):
        cds_list = tx["cds_coords"]
        if not cds_list:
            continue

        chrom = tx["chrom"] or info["chrom"]
        strand = tx["strand"] or info["strand"]
        cds_list = sorted(cds_list, key=lambda x: x[0])
        spliced_dna = ""
        for c_start, c_end in cds_list:
            spliced_dna += extract_sequence_chunk(fasta_file, chrom, c_start, c_end)

        if strand == "-":
            spliced_dna = reverse_complement(spliced_dna)

        protein_seq = translate_dna_to_protein(spliced_dna)
        if protein_seq.endswith("*"):
            protein_seq = protein_seq[:-1]

        records.append({
            "gene_id": gene_id,
            "tx_id": tx_id,
            "chrom": chrom,
            "strand": strand,
            "exon_count": len(cds_list),
            "cds_length": len(spliced_dna),
            "protein": protein_seq,
        })
    return records


def write_protein_records(handle, records):
    for record in records:
        header = (
            f">{record['gene_id']}|{record['tx_id']}_Protein "
            f"{record['chrom']} strand:{record['strand']} "
            f"cds_parts:{record['exon_count']} cds_length:{record['cds_length']}bp "
            f"length:{len(record['protein'])}aa"
        )
        handle.write(format_fasta(header, record["protein"]))


def process_and_export(fasta_file, found_genes, output_mode, protein_mode, base_dir):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    out_dna = os.path.join(base_dir, f"Gene_DNA_{timestamp}.fasta")
    out_prot_all = os.path.join(base_dir, f"CDS_Protein_all_transcripts_{timestamp}.fasta")
    out_prot_longest = os.path.join(base_dir, f"CDS_Protein_longest_transcript_{timestamp}.fasta")

    f_dna = open(out_dna, "w", encoding="utf-8") if output_mode in {"1", "3"} else None
    f_prot_all = open(out_prot_all, "w", encoding="utf-8") if protein_mode in {"1", "3"} else None
    f_prot_longest = open(out_prot_longest, "w", encoding="utf-8") if protein_mode in {"2", "3"} else None

    print("\n正在从参考基因组中提取并转换序列...")

    try:
        for gene_id, info in sorted(found_genes.items()):
            chrom = info["chrom"]
            strand = info["strand"]

            if f_dna:
                if not info.get("gene_coord"):
                    print(f"警告: {gene_id} 缺少基因坐标，跳过 DNA 提取。")
                else:
                    g_start, g_end = info["gene_coord"]
                    raw_dna = extract_sequence_chunk(fasta_file, chrom, g_start, g_end)
                    if raw_dna:
                        final_dna = reverse_complement(raw_dna) if strand == "-" else raw_dna
                        direction = "minus_strand_RC" if strand == "-" else "plus_strand"
                        header = f">{gene_id}_DNA {chrom}:{g_start}-{g_end} strand:{strand} info:{direction}"
                        f_dna.write(format_fasta(header, final_dna))

            if protein_mode:
                protein_records = build_protein_records(fasta_file, gene_id, info)
                if not protein_records:
                    print(f"警告: {gene_id} 没有可用 CDS 注释，已跳过蛋白提取。")
                    continue

                if f_prot_all:
                    write_protein_records(f_prot_all, protein_records)
                if f_prot_longest:
                    longest = max(
                        protein_records,
                        key=lambda r: (len(r["protein"].replace("*", "")), r["cds_length"], r["tx_id"]),
                    )
                    write_protein_records(f_prot_longest, [longest])
    finally:
        if f_dna:
            f_dna.close()
        if f_prot_all:
            f_prot_all.close()
        if f_prot_longest:
            f_prot_longest.close()

    print("\n提取与翻译完成。")
    if output_mode in {"1", "3"}:
        print(f"  基因核苷酸文件: {out_dna}")
    if protein_mode in {"1", "3"}:
        print(f"  所有转录本蛋白文件: {out_prot_all}")
    if protein_mode in {"2", "3"}:
        print(f"  最长转录本蛋白文件: {out_prot_longest}")


def run_pipeline():
    base_dir = get_base_dir()

    fasta_list = (
        glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True)
        + glob.glob(os.path.join(base_dir, "**", "*.fa"), recursive=True)
        + glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
    )
    selected_fasta = interactive_select(fasta_list, "参考基因组 (.fasta/.fa/.fna)", base_dir)
    if not selected_fasta:
        return
    ref_fasta = selected_fasta[0]

    if not os.path.exists(ref_fasta + ".fai"):
        print("\n正在建立参考基因组索引...")
        subprocess.run(["samtools", "faidx", ref_fasta], check=False)

    anno_list = (
        glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)
        + glob.glob(os.path.join(base_dir, "**", "*.gtf"), recursive=True)
    )
    selected_anno = interactive_select(anno_list, "基因注释文件 (.gff/.gff3/.gtf)", base_dir)
    if not selected_anno:
        return
    gff_file = selected_anno[0]

    target_genes, output_mode, protein_mode = get_target_genes(base_dir)

    found_genes = parse_gff_for_targets(gff_file, target_genes)
    if not found_genes:
        print("没有找到任何有效基因记录，程序终止。")
        return

    process_and_export(ref_fasta, found_genes, output_mode, protein_mode, base_dir)


if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n用户中断，程序退出。")
