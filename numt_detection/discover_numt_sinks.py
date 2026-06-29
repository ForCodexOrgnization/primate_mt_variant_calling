#!/usr/bin/env python3
import argparse
import pysam
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Discover candidate NUMT sink loci from mt-related reads remapped to nuclear genome.")
    p.add_argument("--bam", required=True, help="BAM of mt-related reads remapped to nuclear-only reference")
    p.add_argument("--sample", required=True, help="Sample name")
    p.add_argument("--min-mapq", type=int, default=20, help="Minimum read MAPQ")
    p.add_argument("--min-depth", type=int, default=3, help="Minimum pileup depth to consider a position active")
    p.add_argument("--min-reads", type=int, default=5, help="Minimum unique reads supporting a final interval")
    p.add_argument("--min-len", type=int, default=100, help="Minimum interval length")
    p.add_argument("--merge-gap", type=int, default=50, help="Merge nearby active blocks within this gap")
    p.add_argument("--pad", type=int, default=500, help="Pad final intervals by this many bases on each side")
    p.add_argument("--bed-out", required=True)
    p.add_argument("--tsv-out", required=True)
    return p.parse_args()

def get_active_blocks(bam, min_mapq, min_depth):
    """
    Build active blocks from per-base depth >= min_depth using a sweep-line.

    This avoids a whole-reference pileup walk, which can be very slow on
    fragmented or very large references.
    """
    # chrom -> position event delta (+1 on interval start, -1 on end)
    events = defaultdict(lambda: defaultdict(int))

    for aln in bam.fetch(until_eof=False):
        if aln.is_unmapped:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        # aligned reference blocks from CIGAR (includes deletions in span)
        for start0, end0 in aln.get_blocks():
            if end0 <= start0:
                continue
            ev = events[aln.reference_name]
            ev[start0] += 1
            ev[end0] -= 1

    active = defaultdict(list)
    for chrom, ev in events.items():
        depth = 0
        block_start = None

        for pos0 in sorted(ev):
            prev_depth = depth
            depth += ev[pos0]

            # threshold crossing upward starts a block
            if prev_depth < min_depth and depth >= min_depth:
                block_start = pos0
            # threshold crossing downward ends a block
            elif prev_depth >= min_depth and depth < min_depth and block_start is not None:
                active[chrom].append((block_start, pos0))
                block_start = None

        # no trailing close is needed because block ends at an event boundary

    return active

def merge_blocks(blocks, merge_gap):
    """
    Merge nearby blocks.
    """
    merged = []
    if not blocks:
        return merged
    blocks = sorted(blocks)
    cur_s, cur_e = blocks[0]
    for s, e in blocks[1:]:
        if s <= cur_e + merge_gap:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def count_interval_reads(bam, chrom, start0, end0, min_mapq):
    """
    Count unique read names overlapping interval.
    """
    qnames = set()
    mapqs = []
    for aln in bam.fetch(chrom, start0, end0):
        if aln.is_unmapped:
            continue
        if aln.mapping_quality < min_mapq:
            continue
        qnames.add(aln.query_name)
        mapqs.append(aln.mapping_quality)
    nreads = len(qnames)
    mean_mapq = round(sum(mapqs) / len(mapqs), 2) if mapqs else 0.0
    return nreads, mean_mapq

def get_reference_lengths(bam):
    return {ref: bam.get_reference_length(ref) for ref in bam.references}

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    ref_lens = get_reference_lengths(bam)

    active = get_active_blocks(bam, args.min_mapq, args.min_depth)

    rows = []
    bed_rows = []

    for chrom, blocks in active.items():
        merged = merge_blocks(blocks, args.merge_gap)
        chrom_len = ref_lens[chrom]

        for s, e in merged:
            interval_len = e - s
            if interval_len < args.min_len:
                continue

            nreads, mean_mapq = count_interval_reads(bam, chrom, s, e, args.min_mapq)
            if nreads < args.min_reads:
                continue

            pad_s = max(0, s - args.pad)
            pad_e = min(chrom_len, e + args.pad)

            rows.append({
                "sample": args.sample,
                "chrom": chrom,
                "start0": s,
                "end0": e,
                "length": interval_len,
                "nreads": nreads,
                "mean_mapq": mean_mapq,
                "padded_start0": pad_s,
                "padded_end0": pad_e,
                "padded_length": pad_e - pad_s
            })

            bed_rows.append((chrom, pad_s, pad_e, args.sample, nreads, mean_mapq))

    # sort
    rows.sort(key=lambda x: (x["chrom"], x["start0"], x["end0"]))
    bed_rows.sort(key=lambda x: (x[0], x[1], x[2]))

    with open(args.tsv_out, "w") as out:
        out.write("\t".join([
            "sample", "chrom", "start0", "end0", "length",
            "nreads", "mean_mapq",
            "padded_start0", "padded_end0", "padded_length"
        ]) + "\n")
        for r in rows:
            out.write("\t".join(map(str, [
                r["sample"], r["chrom"], r["start0"], r["end0"], r["length"],
                r["nreads"], r["mean_mapq"],
                r["padded_start0"], r["padded_end0"], r["padded_length"]
            ])) + "\n")

    with open(args.bed_out, "w") as out:
        for chrom, s, e, sample, nreads, mean_mapq in bed_rows:
            out.write(f"{chrom}\t{s}\t{e}\t{sample};nreads={nreads};meanMAPQ={mean_mapq}\n")

    bam.close()

if __name__ == "__main__":
    main()