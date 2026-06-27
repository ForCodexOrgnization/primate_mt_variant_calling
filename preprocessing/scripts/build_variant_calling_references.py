#!/usr/bin/env python3
"""Build final variant-calling reference packages after in-house score.

Manual checks for chrM consistency after rebuilding selected rows:

# Show rows where Ref_whole chrM differs from Ref_chrM
awk -F'\t' '
NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
$h["check_status"]!="pass" {
  print $h["target_species"], $h["REF_TYPE"], $h["expected_chrM_mode"], $h["whole_chrM_len"], $h["Ref_chrM_len"], $h["check_message"]
}
' results/preprocessing/reports/variant_reference_integrity_check.tsv

# Compare a single species manually
sp=Aotus_azarae
grep -c "^>chrM$" references/variant_calling/Ref_whole/${sp}.fa
diff -q \
  <(awk '/^>/{p=($0==">chrM"); next} p' references/variant_calling/Ref_whole/${sp}.fa) \
  <(awk '/^>/{p=($0==">chrM"); next} p' references/variant_calling/Ref_chrM/${sp}.fa)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
FASTA_WIDTH = 80
LEFT_MARGIN = 575
RIGHT_MARGIN = 545
CHRM_KEYWORDS = re.compile(r"(^|[^A-Za-z0-9])(chrM|MT)([^A-Za-z0-9]|$)|mitochondr|mitochondria|mtdna", re.I)
FINAL_MASK_PRIORITIES = {
    "C_likely_comp_FINAL_minimal_mask", "C_likely_comp_FINAL_minimal_mask_target_not_reached",
    "C_Ambiguous_FINAL_minimal_mask", "C_Ambiguous_FINAL_minimal_mask_target_not_reached",
    "A_FINAL_minimal_non_chrM_mask", "A_FINAL_minimal_non_chrM_mask_target_not_reached",
}
REQUIRED_REF_COLUMNS = ["target_species", "wg_fasta_path", "wg_fai_path", "chrM_fasta_path", "chrM_fai_path", "chrM_reference_context", "reference_pairing_status", "final_reference_strategy", "final_wg_ref_species", "final_wg_assembly_accession", "final_chrM_species", "final_chrM_accession"]
REQUIRED_SCORE_COLUMNS = ["Species", "REF_TYPE", "MTLIKE_PATTERN", "ValidAnnotatedMitoContig", "HasValidAnnotatedMito", "MaskPriority", "MinimalMaskBED", "FullMaskBED", "CandidateTSV"]
MANIFEST_COLUMNS = ["target_species", "safe_species_id", "REF_TYPE", "MTLIKE_PATTERN", "ValidAnnotatedMitoContig", "HasValidAnnotatedMito", "MaskPriority", "numt_mask_bed", "numt_mask_applied_to_whole_ref", "whole_fasta", "whole_fai", "whole_dict", "chrM_fasta", "chrM_fai", "chrM_dict", "chrM_shift_fasta", "chrM_shift_fai", "chrM_shift_dict", "non_control_interval", "control_region_shifted_interval", "interval_non_control_start", "interval_non_control_end", "interval_control_shifted_start", "interval_control_shifted_end", "shift_back_chain", "wg_fasta_source", "chrM_fasta_source", "chrM_reference_context", "reference_pairing_status", "final_reference_strategy", "final_wg_ref_species", "final_wg_assembly_accession", "final_chrM_species", "final_chrM_accession", "build_status", "build_message"]


def rel(path: Optional[Path]) -> str:
    if not path:
        return ""
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except Exception:
        return str(path)


def resolve_path(value: str) -> Path:
    p = Path((value or "").strip())
    return p if p.is_absolute() else REPO_ROOT / p


def safe_species_id(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", (value or "unknown").strip()).strip("_")
    return safe or "unknown"


def normalize_species_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")


def same_species(a: str, b: str) -> bool:
    norm_a = normalize_species_name(a)
    return norm_a == normalize_species_name(b) and norm_a != ""


def is_missing(value: object) -> bool:
    return str(value or "").strip() in {"", "NA", "N/A", "nan", "None", "null"}


def read_tsv(path: Path, required: Sequence[str]) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing_or_empty_tsv:{rel(path)}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{rel(path)} missing required columns: {','.join(missing)}")
        return [{k: (v or "") for k, v in row.items()} for row in reader]


def fasta_iter(path: Path) -> Iterable[Tuple[str, str, str]]:
    name = desc = None
    seq: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, desc or name, "".join(seq)
                desc = line[1:].strip()
                name = desc.split()[0] if desc else ""
                seq = []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        yield name, desc or name, "".join(seq)


def write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), FASTA_WIDTH):
                out.write(seq[i:i + FASTA_WIDTH] + "\n")


def select_chrm_sequence(chrm_fasta: Path) -> Tuple[str, List[str]]:
    records = list(fasta_iter(chrm_fasta))
    if not records:
        raise ValueError("empty_chrM_fasta")
    if len(records) == 1:
        return records[0][2], []
    candidates = [(name, seq) for name, desc, seq in records if CHRM_KEYWORDS.search(desc) or CHRM_KEYWORDS.search(name)]
    if len(candidates) == 1:
        return candidates[0][1], ["multiple_chrM_records_selected_by_keyword"]
    raise ValueError("ambiguous_multi_record_chrM_fasta")


def bed_intervals(path: Optional[Path]) -> List[Tuple[str, int, int]]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    intervals = []
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((fields[0], start, end))
    return intervals


def mask_sequence(seq: str, intervals: Sequence[Tuple[int, int]]) -> str:
    chars = list(seq)
    for start, end in intervals:
        start, end = max(0, start), min(len(chars), end)
        if end > start:
            chars[start:end] = "N" * (end - start)
    return "".join(chars)


def final_mask_candidate(score: Dict[str, str]) -> Optional[Path]:
    minimal = score.get("MinimalMaskBED", "")
    if not is_missing(minimal):
        min_path = resolve_path(minimal)
        stem = min_path.name
        guesses = []
        if ".minimal_chrMcover.bed" in stem:
            guesses.append(min_path.with_name(stem.replace(".minimal_chrMcover.bed", ".FINAL_numt_mask.bed")))
        guesses.append(min_path.with_name(re.sub(r"\.bed$", ".FINAL_numt_mask.bed", stem)))
        guesses.append(min_path)
        for guess in guesses:
            if guess.exists():
                return guess
        return min_path
    return None


def should_apply_mask(score: Dict[str, str], mask_ref_types: set, bed: Optional[Path], intervals: Sequence[Tuple[str, int, int]]) -> Tuple[bool, str]:
    if score.get("REF_TYPE", "") not in mask_ref_types:
        return False, ""
    if score.get("MaskPriority", "") not in FINAL_MASK_PRIORITIES:
        return False, "mask priority is not final minimal mask"
    if not bed or not bed.exists():
        return False, "missing mask BED"
    if not intervals:
        return False, "mask BED has no intervals"
    return True, ""


def build_whole(
    wg: Path,
    out: Path,
    chrm_seq: str,
    valid_contig: str,
    has_valid_mito: str,
    intervals: Sequence[Tuple[str, int, int]],
    apply_mask: bool,
    target_species: str = "",
    final_chrM_species: str = "",
) -> List[str]:
    warnings: List[str] = []
    by_contig: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for contig, start, end in intervals:
        by_contig[contig].append((start, end))
    embedded = str(has_valid_mito).strip() == "1" and not is_missing(valid_contig)
    selected_chrm_is_target_species = same_species(final_chrM_species, target_species)
    valid_contig = valid_contig.strip()
    chrm_count = 0
    seen = set()
    records = []
    for name, _desc, seq in fasta_iter(wg):
        seen.add(name)
        out_name = "chrM" if embedded and name == valid_contig else name
        if name == "chrM" and out_name != "chrM":
            warnings.append("WG_contig_named_chrM_not_ValidAnnotatedMitoContig")
        if out_name == "chrM":
            chrm_count += 1
            if embedded and name == valid_contig:
                if selected_chrm_is_target_species:
                    if seq != chrm_seq:
                        warnings.append("embedded_chrM_sequence_differed_from_selected_target_chrM_replaced")
                    seq = chrm_seq
                else:
                    warnings.append("selected_chrM_not_target_species_embedded_chrM_not_replaced")
        if apply_mask and out_name != "chrM" and name in by_contig:
            seq = mask_sequence(seq, by_contig[name])
        records.append((out_name, seq))
    if apply_mask:
        missing = sorted(set(by_contig) - seen)
        if missing:
            warnings.append("mask BED contig not found:" + ",".join(missing[:10]))
    if not embedded:
        if any(name == "chrM" for name, _ in records):
            warnings.append("existing chrM contig present before append")
        if not selected_chrm_is_target_species:
            warnings.append("selected_chrM_not_target_species_chrM_appended")
        records.append(("chrM", chrm_seq))
        chrm_count += 1
    write_fasta(out, records)
    if chrm_count != 1:
        raise ValueError(f"Ref_whole_chrM_count_not_one:{chrm_count}")
    return warnings


def fasta_lengths(path: Path) -> Dict[str, int]:
    return {name: len(seq) for name, _desc, seq in fasta_iter(path)}


def md5_sequence(seq: str) -> str:
    return hashlib.md5(seq.encode("ascii")).hexdigest()


def write_intervals(sid: str, chrm_fa: Path, shift_fa: Path, shift: int, interval_dir: Path) -> Tuple[Path, Path, int, int, int, int]:
    seq = next(fasta_iter(chrm_fa))[2]
    shift_seq = next(fasta_iter(shift_fa))[2]
    L = len(seq)
    m5_by_fasta = {chrm_fa: md5_sequence(seq), shift_fa: md5_sequence(shift_seq)}
    nc_start = LEFT_MARGIN + 1
    nc_end = L - RIGHT_MARGIN

    def shifted_pos(x: int) -> int:
        return ((x - 1 - shift) % L + L) % L + 1

    ctrl_start, ctrl_end = sorted((shifted_pos(nc_end + 1), shifted_pos(nc_start - 1)))
    non = interval_dir / f"{sid}_non_control_region.interval_list"
    ctrl = interval_dir / f"{sid}_control_region_shifted.interval_list"
    for path, fasta, start, end in [(non, chrm_fa, nc_start, nc_end), (ctrl, shift_fa, ctrl_start, ctrl_end)]:
        path.write_text(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrM\tLN:{L}\tM5:{m5_by_fasta[fasta]}\tUR:file://{fasta.resolve()}\nchrM\t{start}\t{end}\t+\t.\n")
    print(f"Intervals for {sid}: non_control={nc_start}-{nc_end}; control_shifted={ctrl_start}-{ctrl_end}", file=sys.stderr)
    return non, ctrl, nc_start, nc_end, ctrl_start, ctrl_end


def write_chain(path: Path, L: int, shift: int) -> None:
    S = shift % L
    left = L - S
    path.write_text(f"chain {left} chrM {L} + 0 {left} chrM {L} + {S} {L} 1\n{left}\n\nchain {S} chrM {L} + {left} {L} chrM {L} + 0 {S} 2\n{S}\n\n")


def expected_indexes(fasta: Path) -> List[Path]:
    return [Path(str(fasta) + ext) for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]] + [fasta.with_suffix(".dict")]


def index_fasta(fasta: Path, args: argparse.Namespace) -> None:
    expected = expected_indexes(fasta)
    have = [p.exists() for p in expected]
    if all(have) and not args.force:
        return
    for p in expected:
        if p.exists():
            p.unlink()
    subprocess.run([args.samtools, "faidx", str(fasta)], check=True)
    subprocess.run([args.bwa, "index", str(fasta)], check=True)
    subprocess.run([args.gatk, "CreateSequenceDictionary", "-R", str(fasta), "-O", str(fasta.with_suffix(".dict"))], check=True)


def manifest_paths(base: Dict[str, str], whole: Path, chrm: Path, shift: Path, non: Path, ctrl: Path, chain: Path) -> None:
    base.update({"whole_fasta": rel(whole), "whole_fai": rel(Path(str(whole) + ".fai")), "whole_dict": rel(whole.with_suffix(".dict")), "chrM_fasta": rel(chrm), "chrM_fai": rel(Path(str(chrm) + ".fai")), "chrM_dict": rel(chrm.with_suffix(".dict")), "chrM_shift_fasta": rel(shift), "chrM_shift_fai": rel(Path(str(shift) + ".fai")), "chrM_shift_dict": rel(shift.with_suffix(".dict")), "non_control_interval": rel(non), "control_region_shifted_interval": rel(ctrl), "shift_back_chain": rel(chain)})


def expected_package_paths(sid: str, dirs: Dict[str, Path]) -> Tuple[Path, Path, Path, Path, Path, Path]:
    return (dirs["Ref_whole"] / f"{sid}.fa", dirs["Ref_chrM"] / f"{sid}.fa", dirs["Ref_chrM_shift"] / f"{sid}.fa", dirs["interval"] / f"{sid}_non_control_region.interval_list", dirs["interval"] / f"{sid}_control_region_shifted.interval_list", dirs["shift_back_chain"] / f"{sid}_ShiftBack.chain")


def package_complete(paths: Sequence[Path]) -> bool:
    whole, chrm, shift, non, ctrl, chain = paths
    required = [whole, chrm, shift, non, ctrl, chain]
    for fasta in [whole, chrm, shift]:
        required.extend(expected_indexes(fasta))
    return all(p.exists() and p.stat().st_size > 0 for p in required)


def successful_manifest_rows(manifest: Path) -> Dict[str, Dict[str, str]]:
    if not manifest.exists() or manifest.stat().st_size == 0:
        return {}
    try:
        rows = read_tsv(manifest, MANIFEST_COLUMNS)
    except Exception as exc:
        print(f"WARNING: could not read existing manifest for skip checks: {exc}", file=sys.stderr)
        return {}
    return {r["safe_species_id"]: r for r in rows if r.get("build_status") == "success" and r.get("safe_species_id")}


def validate_manifest(manifest: Path) -> int:
    rows = read_tsv(manifest, MANIFEST_COLUMNS)
    success = [r for r in rows if r.get("build_status") == "success"]
    failed = [r for r in rows if r.get("build_status") != "success"]
    errors = []
    for r in success:
        sid = r["safe_species_id"]
        whole, chrm, shift = map(resolve_path, [r["whole_fasta"], r["chrM_fasta"], r["chrM_shift_fasta"]])
        wlens, clens, slens = fasta_lengths(whole), fasta_lengths(chrm), fasta_lengths(shift)
        if list(clens) != ["chrM"] or list(slens) != ["chrM"] or len(clens) != 1 or len(slens) != 1:
            errors.append(f"{sid}: Ref_chrM/shift record name/count invalid")
        if clens.get("chrM") != slens.get("chrM"):
            errors.append(f"{sid}: chrM and shifted lengths differ")
        if sum(1 for n in wlens if n == "chrM") != 1:
            errors.append(f"{sid}: Ref_whole does not contain exactly one chrM")
        for fasta in [whole, chrm, shift]:
            missing = [rel(p) for p in expected_indexes(fasta) if not p.exists()]
            if missing:
                errors.append(f"{sid}: missing indexes for {rel(fasta)}: {','.join(missing)}")
    print(f"Validation summary: successful={len(success)} failed_or_skipped={len(failed)}")
    table = Counter((r.get("REF_TYPE", ""), r.get("MaskPriority", ""), r.get("numt_mask_applied_to_whole_ref", "")) for r in rows)
    if table:
        print("REF_TYPE\tMaskPriority\tnumt_mask_applied_to_whole_ref\tcount")
        for (ref_type, priority, applied), count in sorted(table.items()):
            print(f"{ref_type}\t{priority}\t{applied}\t{count}")
    for error in errors:
        print("VALIDATION_ERROR\t" + error, file=sys.stderr)
    return 1 if errors else 0


def build_one_reference(ref: Dict[str, str], sid: str, score: Dict[str, str], dirs: Dict[str, Path], existing_success: Dict[str, Dict[str, str]], mask_ref_types: set, args: argparse.Namespace) -> Dict[str, str]:
    sp = ref["target_species"]
    base = {c: "" for c in MANIFEST_COLUMNS}
    base.update({"target_species": sp, "safe_species_id": sid, "wg_fasta_source": ref["wg_fasta_path"], "chrM_fasta_source": ref["chrM_fasta_path"], "REF_TYPE": score.get("REF_TYPE", ""), "MTLIKE_PATTERN": score.get("MTLIKE_PATTERN", ""), "ValidAnnotatedMitoContig": score.get("ValidAnnotatedMitoContig", ""), "HasValidAnnotatedMito": score.get("HasValidAnnotatedMito", ""), "MaskPriority": score.get("MaskPriority", ""), "chrM_reference_context": ref.get("chrM_reference_context", ""), "reference_pairing_status": ref.get("reference_pairing_status", ""), "final_reference_strategy": ref.get("final_reference_strategy", ""), "final_wg_ref_species": ref.get("final_wg_ref_species", ""), "final_wg_assembly_accession": ref.get("final_wg_assembly_accession", ""), "final_chrM_species": ref.get("final_chrM_species", ""), "final_chrM_accession": ref.get("final_chrM_accession", "")})
    messages: List[str] = []
    try:
        whole_fa, chrm_fa, shift_fa, non, ctrl, chain = expected_package_paths(sid, dirs)
        if sid in existing_success and package_complete([whole_fa, chrm_fa, shift_fa, non, ctrl, chain]) and not args.force:
            base.update(existing_success[sid])
            base.update({"build_status": "success", "build_message": "skipped_existing_manifest_success_and_outputs"})
            manifest_paths(base, whole_fa, chrm_fa, shift_fa, non, ctrl, chain)
            return base
        wg, chrm_src = resolve_path(ref["wg_fasta_path"]), resolve_path(ref["chrM_fasta_path"])
        if not chrm_src.exists():
            raise FileNotFoundError("missing_chrM_fasta")
        if not wg.exists():
            raise FileNotFoundError("missing_wg_fasta")
        mask_bed = final_mask_candidate(score)
        intervals = bed_intervals(mask_bed)
        apply_mask, mask_msg = should_apply_mask(score, mask_ref_types, mask_bed, intervals)
        chrm_seq, chrm_warn = select_chrm_sequence(chrm_src)
        messages.extend(chrm_warn)
        write_fasta(chrm_fa, [("chrM", chrm_seq)])
        S = args.shift % len(chrm_seq)
        write_fasta(shift_fa, [("chrM", chrm_seq[S:] + chrm_seq[:S])])
        if mask_msg and score.get("REF_TYPE", "") in mask_ref_types:
            messages.append(mask_msg)
        messages.extend(build_whole(wg=wg, out=whole_fa, chrm_seq=chrm_seq, valid_contig=score.get("ValidAnnotatedMitoContig", ""), has_valid_mito=score.get("HasValidAnnotatedMito", ""), intervals=intervals, apply_mask=apply_mask, target_species=ref.get("target_species", ""), final_chrM_species=ref.get("final_chrM_species", "")))
        non, ctrl, nc_start, nc_end, ctrl_start, ctrl_end = write_intervals(sid, chrm_fa, shift_fa, args.shift, dirs["interval"])
        base.update({"interval_non_control_start": str(nc_start), "interval_non_control_end": str(nc_end), "interval_control_shifted_start": str(ctrl_start), "interval_control_shifted_end": str(ctrl_end)})
        write_chain(chain, len(chrm_seq), args.shift)
        for fasta in [whole_fa, chrm_fa, shift_fa]:
            index_fasta(fasta, args)
        base.update({"numt_mask_bed": rel(mask_bed), "numt_mask_applied_to_whole_ref": "yes" if apply_mask else "no", "build_status": "success", "build_message": ";".join(messages) or "ok"})
        manifest_paths(base, whole_fa, chrm_fa, shift_fa, non, ctrl, chain)
    except Exception as exc:
        base.update({"build_status": "failed", "build_message": ";".join(messages + [str(exc)])})
    return base


def build(args: argparse.Namespace) -> Path:
    out_root = resolve_path(args.out_root)
    dirs = {name: out_root / name for name in ["Ref_whole", "Ref_chrM", "Ref_chrM_shift", "interval", "shift_back_chain"]}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    existing_manifest = out_root / "variant_calling_reference_manifest.tsv"
    existing_success = successful_manifest_rows(existing_manifest) if not args.force else {}
    refs = read_tsv(resolve_path(args.ref_inputs), REQUIRED_REF_COLUMNS)
    scores = read_tsv(resolve_path(args.score), REQUIRED_SCORE_COLUMNS)
    score_by_species = {r["Species"]: r for r in scores}
    seen, unique = set(), []
    for r in refs:
        key = (r["target_species"], r["wg_fasta_path"], r["chrM_fasta_path"])
        if key not in seen:
            seen.add(key)
            unique.append(r)
    species_counts = Counter(r["target_species"] for r in unique)
    species_ord = defaultdict(int)
    jobs = []
    mask_ref_types = {x.strip() for x in args.mask_ref_types.split(",") if x.strip()}
    for ref in unique:
        sp = ref["target_species"]
        species_ord[sp] += 1
        sid = safe_species_id(sp)
        if species_counts[sp] > 1:
            sid = f"{sid}_ref{species_ord[sp]:04d}"
        jobs.append((ref, sid, score_by_species.get(sp, {})))
    threads = max(1, args.threads)
    if args.job_index:
        if args.job_index < 1 or args.job_index > len(jobs):
            raise SystemExit(f"--job-index {args.job_index} outside available reference package range 1..{len(jobs)}")
        rows = [build_one_reference(*jobs[args.job_index - 1], dirs, existing_success, mask_ref_types, args)]
    elif threads == 1 or len(jobs) <= 1:
        rows = [build_one_reference(ref, sid, score, dirs, existing_success, mask_ref_types, args) for ref, sid, score in jobs]
    else:
        print(f"Building {len(jobs)} variant-reference packages with {threads} workers", file=sys.stderr)
        rows = [None] * len(jobs)
        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_index = {executor.submit(build_one_reference, ref, sid, score, dirs, existing_success, mask_ref_types, args): i for i, (ref, sid, score) in enumerate(jobs)}
            for future in as_completed(future_to_index):
                rows[future_to_index[future]] = future.result()
    manifest = resolve_path(args.manifest_output) if args.manifest_output else out_root / "variant_calling_reference_manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {rel(manifest)}")
    validate_manifest(manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-inputs", default="references/manifests/in_house_score_reference_inputs.tsv")
    parser.add_argument("--score", "--score-file", dest="score", default="results/preprocessing/in_house_score/merged_in_house_score.tsv")
    parser.add_argument("--out-root", "--outdir", dest="out_root", default="references/variant_calling")
    parser.add_argument("--shift", type=int, default=8000)
    parser.add_argument("--mask-ref-types", default="#C-likely_comp,#C-Ambiguous,#A")
    parser.add_argument("--threads", type=int, default=int(os.environ.get("VARIANT_REFERENCE_THREADS", "1")), help="Number of reference packages to build in parallel")
    parser.add_argument("--job-index", type=int, default=0, help="Build only the 1-based reference package index, for Slurm array tasks")
    parser.add_argument("--manifest-output", default="", help="Write manifest rows to this TSV instead of the default final manifest")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--samtools", default=os.environ.get("SAMTOOLS_COMMAND", "samtools"))
    parser.add_argument("--bwa", default=os.environ.get("BWA_COMMAND", "bwa"))
    parser.add_argument("--gatk", default=os.environ.get("GATK_COMMAND", "gatk"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = resolve_path(args.out_root) / "variant_calling_reference_manifest.tsv"
    if args.validate_only:
        return validate_manifest(manifest)
    build(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
