#!/usr/bin/env python3
"""Prepare full nuclear regions and summarize mosdepth output for round2 mtCN."""

import argparse
import csv
import gzip
import math
from pathlib import Path
from statistics import median


def read_fai(path, mt_contig):
    with open(path) as handle:
        return [(f[0], int(f[1])) for line in handle if len(f := line.rstrip().split("\t")) >= 2 and f[0] != mt_contig]


def choose_regions(contigs, initial, target):
    total = sum(length for _, length in contigs)
    thresholds = [initial] + [x for x in (50000, 20000, 10000, 5000, 1000, 0) if x < initial]
    for threshold in dict.fromkeys(thresholds):
        selected = [(name, length) for name, length in contigs if length >= threshold]
        fraction = sum(length for _, length in selected) / total if total else 0.0
        if fraction >= target or threshold == 0:
            return threshold, selected, total, fraction
    raise AssertionError("threshold list did not include zero")


def read_cram_header(path):
    contigs = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith("@SQ\t"):
                continue
            tags = dict(field.split(":", 1) for field in line.rstrip().split("\t")[1:] if ":" in field)
            if "SN" in tags and "LN" in tags:
                contigs[tags["SN"]] = (int(tags["LN"]), tags.get("M5", ""))
    return contigs


def prepare(args):
    contigs = read_fai(args.fai, args.mt_contig)
    threshold, selected, total, fraction = choose_regions(contigs, args.min_length, args.target_fraction)
    if not selected:
        raise SystemExit("ERROR: no non-mt nuclear contigs found in reference")
    Path(args.bed).write_text("".join(f"{name}\t0\t{length}\n" for name, length in selected))
    cram = read_cram_header(args.cram_header)
    failures = []
    with open(args.validation, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["contig", "reference_length", "cram_length", "cram_M5", "status"])
        for name, length in selected:
            if name not in cram:
                row = [name, length, "", "", "MISSING_IN_CRAM"]
            elif cram[name][0] != length:
                row = [name, length, cram[name][0], cram[name][1], "LENGTH_MISMATCH"]
            else:
                row = [name, length, cram[name][0], cram[name][1], "MATCH"]
            out.writerow(row)
            if row[-1] != "MATCH":
                failures.append(row)
    selected_bp = sum(length for _, length in selected)
    with open(args.qc, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample", "ref_name", "nuclear_min_contig_len_used", "nuclear_contigs", "selected_nuclear_bp", "total_non_mt_bp", "selected_bp_fraction"])
        out.writerow([args.sample, args.ref_name, threshold, len(selected), selected_bp, total, f"{fraction:.9f}"])
    if failures:
        details = ", ".join(f"{row[0]}:{row[-1]}" for row in failures[:10])
        raise SystemExit(f"ERROR: nuclear reference is incompatible with CRAM header ({details}); see {args.validation}")


def survival_distribution(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3 and fields[0] == "total":
                values[int(fields[1])] = float(fields[2])
    if not values:
        raise SystemExit(f"ERROR: no total cumulative distribution rows in {path}")
    return values


def survival_at(values, depth):
    return values.get(depth, 0.0 if depth > max(values) else max(v for d, v in values.items() if d >= depth))


def quantile(values, survival_cutoff):
    eligible = [depth for depth, fraction in values.items() if fraction >= survival_cutoff]
    return max(eligible) if eligible else 0


def summarize(args):
    nuc_bases = 0
    weighted_depth = 0.0
    with gzip.open(args.regions, "rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 4 or line.startswith("#"):
                continue
            length = int(fields[2]) - int(fields[1])
            nuc_bases += length
            weighted_depth += float(fields[-1]) * length
    if nuc_bases <= 0:
        raise SystemExit("ERROR: no nuclear coverage rows found in mosdepth output")
    mt_depths = [float(line.rstrip().split("\t")[2]) for line in open(args.mt_depth) if line.strip()]
    if not mt_depths:
        raise SystemExit("ERROR: no mt positions found in round2 standard chrM BAM depth")
    dist = survival_distribution(args.distribution)
    nuc_mean = weighted_depth / nuc_bases
    nuc_median = quantile(dist, .50)
    mt_mean = sum(mt_depths) / len(mt_depths)
    mt_med = median(mt_depths)
    mtcn_mean = 2 * mt_mean / nuc_mean if nuc_mean > 0 else math.nan
    mtcn_med = 2 * mt_med / nuc_median if nuc_median > 0 else math.nan
    qc = next(csv.DictReader(open(args.qc), delimiter="\t"))
    extra_names = ["nuclear_q25", "nuclear_q75", "nuclear_zero_coverage_fraction", "nuclear_fraction_ge_10x", "nuclear_fraction_ge_20x", "nuclear_fraction_ge_30x", "nuclear_selected_contigs", "nuclear_selected_bp", "nuclear_total_non_mt_bp", "nuclear_selected_bp_fraction", "nuclear_min_contig_len_used"]
    header = "sample species ref_name mt_contig mt_coverage_source nuclear_coverage_source mean_mt_coverage mean_nuclear_coverage mean_mtCN mt_mean_coverage nuclear_mean_coverage mtcn_mean mean_formula mt_median_coverage nuclear_median_coverage mtcn_median median_formula".split() + extra_names
    extras = [quantile(dist, .75), quantile(dist, .25), 1-survival_at(dist, 1), survival_at(dist, 10), survival_at(dist, 20), survival_at(dist, 30), qc["nuclear_contigs"], qc["selected_nuclear_bp"], qc["total_non_mt_bp"], qc["selected_bp_fraction"], qc["nuclear_min_contig_len_used"]]
    row = [args.sample, args.species, args.ref_name, args.mt_contig, "round2_standard_chrM_assigned_bam", "wgs_cram_mosdepth_full_nuclear_regions", mt_mean, nuc_mean, mtcn_mean, mt_mean, nuc_mean, mtcn_mean, "2*mt_mean_coverage/nuclear_mean_coverage", mt_med, nuc_median, mtcn_med, "2*mt_median_coverage/nuclear_median_coverage"] + extras
    with open(args.output, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(header)
        out.writerow([f"{x:.6f}" if isinstance(x, float) else x for x in row])


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser("prepare")
    for flag in ("sample", "ref-name", "fai", "mt-contig", "cram-header", "bed", "validation", "qc"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--min-length", type=int, default=50000); p.add_argument("--target-fraction", type=float, default=.90); p.set_defaults(func=prepare)
    p = sub.add_parser("summarize")
    for flag in ("sample", "species", "ref-name", "mt-contig", "regions", "distribution", "mt-depth", "qc", "output"):
        p.add_argument("--" + flag, required=True)
    p.set_defaults(func=summarize)
    args = parser.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
