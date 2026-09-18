#!/usr/bin/env python3
"""Batch extract the longest spliced CDS per gene. Python >= 3.9, no dependencies.

Place this script beside Genome/ and GFF3/, then run python extract_longest_cds.py.
All annotated CDS bases are retained, including partial terminal codons. Internal
CDS phase bases must not be removed: they complete codons across splice junctions.
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from urllib.parse import unquote


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8-sig') if path.suffix == '.gz' else path.open(encoding='utf-8-sig')


def read_genome(path):
    genome = {}
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    genome[name] = ''.join(chunks).upper()
                name = line[1:].split()[0]
                if name in genome:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
                chunks = []
            else:
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                chunks.append(''.join(line.split()))
    if name is not None:
        genome[name] = ''.join(chunks).upper()
    return genome


def read_gff(path):
    parents, genes, cds = {}, set(), defaultdict(set)
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GFF line {number}: expected 9 tab-separated columns')
            seqid, _, kind, start, end, _, strand, phase, raw = fields
            attrs = dict(item.split('=', 1) for item in raw.split(';') if '=' in item)
            identifier = unquote(attrs.get('ID', ''))
            parent_ids = tuple(unquote(p) for p in attrs.get('Parent', '').split(',') if p)
            if identifier:
                parents[identifier] = parent_ids
                if kind in ('gene', 'pseudogene', 'SO:0000704'):
                    genes.add(identifier)
            if kind not in ('CDS', 'SO:0000316'):
                continue
            start, end = int(start), int(end)
            if start < 1 or end < start or strand not in ('+', '-') or phase not in ('0', '1', '2', '.'):
                raise ValueError(f'GFF line {number}: invalid CDS coordinates, strand or phase')
            if not parent_ids:
                raise ValueError(f'GFF line {number}: CDS has no Parent; cannot safely group segments')
            for parent in parent_ids:
                cds[parent].add((unquote(seqid), start, end, strand))
    return parents, genes, cds


def gene_ids(identifier, parents, genes, trail=()):
    if identifier in trail:
        raise ValueError(f'Cycle in GFF Parent hierarchy: {identifier}')
    if identifier in genes:
        return {identifier}
    result = set()
    for parent in parents.get(identifier, ()):
        result.update(gene_ids(parent, parents, genes, trail + (identifier,)))
    return result


COMPLEMENT = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')


def process_species(task):
    species, genome_path, gff_path, output, header_format = task
    try:
        parents, genes, transcripts = read_gff(gff_path)
        if not transcripts:
            raise ValueError('No CDS features found')
        selected = {}
        fallback = 0
        for transcript, segments in sorted(transcripts.items()):
            segments = sorted(segments)
            if len({(s[0], s[3]) for s in segments}) != 1:
                raise ValueError(f'{transcript}: CDS spans multiple contigs or strands')
            if any(a[2] >= b[1] for a, b in zip(segments, segments[1:])):
                raise ValueError(f'{transcript}: overlapping CDS segments')
            length = sum(end - start + 1 for _, start, end, _ in segments)
            targets = gene_ids(transcript, parents, genes)
            if not targets:
                targets = {transcript}
                fallback += 1
            for gene in sorted(targets):
                if gene not in selected or length > selected[gene][0]:
                    selected[gene] = (length, transcript, segments)
        genome = read_genome(genome_path)
        # Validate even unselected transcripts, so annotation mismatches cannot go unnoticed.
        for transcript, segments in transcripts.items():
            for contig, start, end, strand in segments:
                if contig not in genome or end > len(genome[contig]):
                    raise ValueError(f'{transcript}: contig missing or coordinates exceed genome ({contig}:{start}-{end})')
        prefix = species.zfill(6) if species.isdecimal() else species
        fasta_path = output / f'{species}.cds.fna'
        map_path = output / f'{species}.cds.tsv'
        if fasta_path.exists() or map_path.exists():
            raise FileExistsError(f'Output already exists for {species}; choose a new --output directory')
        fasta_tmp = fasta_path.with_suffix('.fna.tmp')
        map_tmp = map_path.with_suffix('.tsv.tmp')
        try:
            with fasta_tmp.open('w', encoding='utf-8', newline='\n') as fasta, map_tmp.open('w', encoding='utf-8', newline='') as mapping:
                writer = csv.writer(mapping, delimiter='\t')
                writer.writerow(['sequence_id', 'species_id', 'gene_id', 'transcript_id', 'cds_length', 'contig', 'strand'])
                for index, (gene, (length, transcript, segments)) in enumerate(sorted(selected.items()), 1):
                    sequence = ''.join(genome[c][start-1:end] for c, start, end, _ in segments)
                    if segments[0][3] == '-':
                        sequence = sequence.translate(COMPLEMENT)[::-1]
                    identifier = f'{prefix}|CDS{index}' if header_format == 'species-cds' else f'{index:06d}'
                    fasta.write(f'>{identifier}\n')
                    for offset in range(0, len(sequence), 80):
                        fasta.write(sequence[offset:offset+80] + '\n')
                    writer.writerow([identifier, species, gene, transcript, length, segments[0][0], segments[0][3]])
            fasta_tmp.replace(fasta_path)
            map_tmp.replace(map_path)
        finally:
            fasta_tmp.unlink(missing_ok=True)
            map_tmp.unlink(missing_ok=True)
        return species, len(selected), fallback, ''
    except Exception as exc:
        return species, 0, 0, str(exc)


def discover(folder, extensions):
    result = {}
    for path in sorted(folder.iterdir()):
        if not path.is_file():
            continue
        name = path.name[:-3] if path.name.endswith('.gz') else path.name
        for extension in extensions:
            if name.lower().endswith(extension):
                species = name[:-len(extension)]
                if species in result:
                    raise ValueError(f'Duplicate species ID in {folder}: {species}')
                result[species] = path
                break
    return result


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', type=Path, default=base / 'Genome')
    parser.add_argument('--gff3', type=Path, default=base / 'GFF3')
    parser.add_argument('--output', type=Path, default=base / 'CDS')
    parser.add_argument('--workers', type=int, default=1, help='Parallel species; each worker loads one genome into RAM')
    parser.add_argument('--header-format', choices=['species-cds', 'six-digit'], default='species-cds', help='six-digit uses >000001 etc.; numbering restarts in each species file')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be >= 1')
    try:
        genomes = discover(args.genome, ('.fna', '.fa', '.fasta'))
        annotations = discover(args.gff3, ('.gff3', '.gff'))
        species_ids = sorted(genomes.keys() & annotations.keys())
        missing = genomes.keys() ^ annotations.keys()
        for species in sorted(missing):
            print(f'[UNPAIRED] {species}', file=sys.stderr)
        if not species_ids:
            raise ValueError('No matching Genome/GFF3 pairs found')
        args.output.mkdir(parents=True, exist_ok=True)
        tasks = [(s, genomes[s], annotations[s], args.output, args.header_format) for s in species_ids]
        failed = bool(missing)
        def report(results):
            nonlocal failed
            for species, count, fallback, error in results:
                if error:
                    failed = True
                    print(f'[FAILED] {species}: {error}', file=sys.stderr, flush=True)
                else:
                    print(f'[OK] {species}: {count} CDS; {fallback} transcripts without gene ancestry', flush=True)
        if args.workers == 1:
            report(map(process_species, tasks))
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                report(executor.map(process_species, tasks))
        return 1 if failed else 0
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
