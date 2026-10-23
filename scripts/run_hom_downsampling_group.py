#!/usr/bin/env python3
"""Prepare group inputs for a single multi-sample HOM Nextflow launch.

No variant calling/downsampling is done here. Reuses existing output paths.
Validates CRAM, candidate BAM, reference, sample-specific NUMT BED intervals.
"""
import argparse
import csv
import os
from pathlib import Path

def fail(message):
    raise SystemExit("ERROR: " + message)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", required=True, type=Path)
    p.add_argument("--depth", required=True, type=int, choices=[30, 50, 100])
    p.add_argument("--seed", required=True, type=int, choices=[101, 202])
    p.add_argument("--cdp", required=True, type=int)
    p.add_argument("--wgs-root", required=True, type=Path)
    p.add_argument("--candidate-meta-root", required=True, type=Path)
    p.add_argument("--prod-numt-dir", required=True, type=Path)
    p.add_argument("--ref-whole", required=True, type=Path)
    p.add_argument("--ref-mt", required=True, type=Path)
    args = p.parse_args()

    base = args.base
    tag = f"ds{args.depth}x_seed{args.seed}_cdp{args.cdp}"
    group = base / "variant_calling_pilot_group_inputs" / tag
    cram_dir = group / "cram"
    numt_dir = group / "numt_bed"
    round1 = base / "variant_calling_pilot" / tag / "round1"
    round2 = base / "variant_calling_pilot" / tag / "round2"
    manifest = base / "candidate_downsampling_jobs.tsv"
    pilot = base / "pilot_samples.tsv"

    if not manifest.is_file() or not pilot.is_file():
        fail(f"Missing manifest or pilot_samples.tsv under {base}")
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with pilot.open(newline="") as handle:
        pilot_rows = list(csv.reader(handle, delimiter="\t"))
    # pilot_samples.tsv is a metadata table, NOT a Nextflow sample.tsv.
    # Its third column may be a numeric QC metric (e.g. 99.988), rather
    # than a reference name. Match the existing single-sample launcher:
    # column 1 = sample ID; column 2 = species AND reference name.
    sample_info = {}
    for row in pilot_rows[1:]:
        if len(row) >= 2 and row[0].strip() and row[1].strip():
            sample = row[0].strip()
            species = row[1].strip()
            if sample in sample_info and sample_info[sample] != (species, species):
                fail(f"Conflicting species for sample {sample} in pilot_samples.tsv")
            sample_info[sample] = (species, species)
    subset = [r for r in rows if int(r["target_depth"]) == args.depth and int(r["seed"]) == args.seed]
    if len(subset) != 10 or len({r["sample"] for r in subset}) != 10:
        fail(f"Expected 10 unique samples for {args.depth}x seed{args.seed}; found {len(subset)}")
    for directory in (group, cram_dir, numt_dir, round1, round2):
        directory.mkdir(parents=True, exist_ok=True)

    def link(src, dst):
        if not src.is_file() or src.stat().st_size == 0:
            fail(f"Missing or empty file: {src}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_symlink():
            if dst.resolve() == src.resolve():
                return
            fail(f"Conflicting symlink (not overwriting): {dst}")
        if dst.exists():
            fail(f"Conflicting existing path (not overwriting): {dst}")
        dst.symlink_to(src.resolve())

    def check_and_make_bed(sample, ref_name):
        destination = numt_dir / f"{sample}.highconf_numt.bed"
        production = args.prod_numt_dir / destination.name
        if production.is_file():
            link(production, destination)
            return
        # Reuse the original per-sample fallback BED if already prepared.
        previous = (base / "variant_calling_pilot_inputs" / sample / tag /
                    "numt_bed_from_extraction" / destination.name)
        if previous.is_file() and previous.stat().st_size:
            link(previous, destination)
            return
        extraction = (args.candidate_meta_root / sample / "round_1" /
                      "candidate_reads" / f"{sample}.target_chrM_plus_NUMT.bed")
        fai = args.ref_whole / f"{ref_name}.fa.fai"
        if not extraction.is_file() or not fai.is_file():
            fail(f"Neither NUMT BED nor extraction BED/reference .fai for {sample}")
        sizes = {}
        with fai.open() as stream:
            for line in stream:
                fields = line.split("\t")
                if len(fields) >= 2:
                    sizes[fields[0]] = int(fields[1])
        lines = []
        with extraction.open() as stream:
            for line in stream:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.split()
                if len(f) < 3:
                    fail(f"Invalid extraction BED line for {sample}: {line.strip()}")
                contig = f[0]
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    fail(f"Invalid coordinates in extraction BED for {sample}")
                if contig == "chrM":
                    continue
                if contig not in sizes or start < 0 or end <= start or end > sizes[contig]:
                    fail(f"Out-of-range nuclear NUMT BED interval for {sample}: {line.strip()}")
                lines.append(f"{contig}\t{start}\t{end}\n")
        if not lines:
            fail(f"No nuclear NUMT intervals for {sample}; manual review required")
        if destination.exists() or destination.is_symlink():
            if destination.is_file() and destination.read_text() == "".join(lines):
                return
            fail(f"Refusing to overwrite different NUMT BED: {destination}")
        with destination.open("x") as handle:
            handle.writelines(lines)
        print(f"NUMT BED derived from extraction BED for {sample}: {len(lines)} intervals")

    sample_lines = []
    for r in subset:
        sample = r["sample"]
        if sample not in sample_info:
            fail(f"Sample {sample} missing from pilot_samples.tsv")
        species, ref_name = sample_info[sample]
        bam = base / "candidate_downsampled" / sample / f"{sample}.candidate.ds{args.depth}x.seed{args.seed}.bam"
        if not bam.is_file() or not (Path(str(bam) + ".bai")).is_file():
            fail(f"Missing candidate BAM/BAI for {sample}: {bam}")
        if not (args.ref_whole / f"{ref_name}.fa").is_file() or not (args.ref_mt / f"{ref_name}.fa").is_file():
            fail(f"Missing reference for {sample}: {ref_name}")
        cram = args.wgs_root / sample / f"{sample}.cram"
        link(cram, cram_dir / sample / "alignment" / f"{sample}.cram")
        link(Path(str(cram) + ".crai"), cram_dir / sample / "alignment" / f"{sample}.cram.crai")
        candidate = round1 / sample / "round_1" / "candidate_reads"
        link(bam, candidate / f"{sample}.with_mates.bam")
        link(Path(str(bam) + ".bai"), candidate / f"{sample}.with_mates.bam.bai")
        check_and_make_bed(sample, ref_name)
        sample_lines.append(f"{sample}\t{species}\t{ref_name}\n")

    samples = group / "sample.tsv"
    contents = "".join(sample_lines)
    if samples.is_file() and samples.read_text() != contents:
        fail(f"Existing group sample TSV differs: {samples}")
    if not samples.exists():
        samples.write_text(contents)
    print(f"GROUP_READY={tag} samples={len(sample_lines)} sample_tsv={samples}")
    print(f"CRAM_DIR={cram_dir}")
    print(f"NUMT_DIR={numt_dir}")
    print(f"ROUND1={round1}")
    print(f"ROUND2={round2}")

if __name__ == "__main__":
    main()
