import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "mtcn_nuclear_metrics.py"


class NuclearMetricsTest(unittest.TestCase):
    def prepare(self, fai_rows, target=.9, header_rows=None):
        tmp = tempfile.TemporaryDirectory()
        root = Path(tmp.name)
        (root / "ref.fa.fai").write_text("".join(f"{n}\t{length}\t0\t0\t0\n" for n, length in fai_rows))
        header_rows = header_rows if header_rows is not None else fai_rows
        (root / "header.sam").write_text("".join(f"@SQ\tSN:{n}\tLN:{length}\tM5:test\n" for n, length in header_rows))
        cmd = [sys.executable, str(SCRIPT), "prepare", "--sample", "S1", "--ref-name", "ref", "--fai", str(root / "ref.fa.fai"), "--mt-contig", "chrM", "--cram-header", str(root / "header.sam"), "--bed", str(root / "regions.bed"), "--validation", str(root / "validation.tsv"), "--qc", str(root / "qc.tsv"), "--min-length", "50000", "--target-fraction", str(target)]
        result = subprocess.run(cmd, text=True, capture_output=True)
        return tmp, root, result

    def test_chromosome_assembly_uses_complete_contigs(self):
        tmp, root, result = self.prepare([("chr1", 100000), ("chr2", 80000), ("chrM", 1000)])
        self.addCleanup(tmp.cleanup)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((root / "regions.bed").read_text(), "chr1\t0\t100000\nchr2\t0\t80000\n")
        with open(root / "qc.tsv") as handle:
            qc = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(qc["selected_bp_fraction"], "1.000000000")

    def test_fragmented_assembly_adapts_threshold(self):
        rows = [("chr1", 60000)] + [(f"scaf{i}", 10000) for i in range(6)] + [("chrM", 1000)]
        tmp, root, result = self.prepare(rows)
        self.addCleanup(tmp.cleanup)
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(root / "qc.tsv") as handle:
            qc = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(qc["nuclear_min_contig_len_used"], "10000")
        self.assertEqual(qc["selected_bp_fraction"], "1.000000000")

    def test_mismatch_fails_before_coverage(self):
        tmp, root, result = self.prepare([("chr1", 100000), ("chrM", 1000)], header_rows=[("chr1", 99999), ("chrM", 1000)])
        self.addCleanup(tmp.cleanup)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("LENGTH_MISMATCH", result.stderr)
        self.assertIn("LENGTH_MISMATCH", (root / "validation.tsv").read_text())

    def test_base_distribution_quantiles_and_compatible_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with gzip.open(root / "regions.gz", "wt") as out:
                out.write("chr1\t0\t100\tchr1\t12.0\n")
            (root / "dist.txt").write_text("total\t0\t1.0\ntotal\t1\t0.9\ntotal\t10\t0.8\ntotal\t11\t0.6\ntotal\t12\t0.5\ntotal\t13\t0.2\ntotal\t20\t0.1\ntotal\t30\t0.05\n")
            (root / "mt.tsv").write_text("chrM\t1\t100\nchrM\t2\t200\n")
            (root / "qc.tsv").write_text("sample\tref_name\tnuclear_min_contig_len_used\tnuclear_contigs\tselected_nuclear_bp\ttotal_non_mt_bp\tselected_bp_fraction\nS1\tref\t50000\t1\t100\t100\t1.0\n")
            cmd = [sys.executable, str(SCRIPT), "summarize", "--sample", "S1", "--species", "sp", "--ref-name", "ref", "--mt-contig", "chrM", "--regions", str(root / "regions.gz"), "--distribution", str(root / "dist.txt"), "--mt-depth", str(root / "mt.tsv"), "--qc", str(root / "qc.tsv"), "--output", str(root / "mtcn.tsv"), "--marker", str(root / "method.tsv")]
            subprocess.run(cmd, check=True)
            with open(root / "mtcn.tsv") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["nuclear_median_coverage"], "12")
            for name in ("mean_mt_coverage", "mean_nuclear_coverage", "mean_mtCN", "mt_mean_coverage", "nuclear_mean_coverage", "mtcn_mean", "mt_median_coverage", "nuclear_median_coverage", "mtcn_median"):
                self.assertIn(name, row)
            with open(root / "method.tsv") as handle:
                marker = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(marker, {
                "sample": "S1", "ref_name": "ref", "nuclear_cov_method": "full_nuclear_regions_v1",
                "nuclear_min_contig_len_used": "50000", "selected_bp_fraction": "1.0",
                "nuclear_mean_coverage": "12.000000", "nuclear_median_coverage": "12",
            })

    def test_completion_marker_is_not_created_on_summary_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with gzip.open(root / "regions.gz", "wt") as out:
                out.write("chr1\t0\t100\tchr1\t12.0\n")
            (root / "dist.txt").write_text("chr1\t0\t1.0\n")  # no required total rows
            (root / "mt.tsv").write_text("chrM\t1\t100\n")
            (root / "qc.tsv").write_text("sample\tref_name\tnuclear_min_contig_len_used\tnuclear_contigs\tselected_nuclear_bp\ttotal_non_mt_bp\tselected_bp_fraction\nS1\tref\t50000\t1\t100\t100\t1.0\n")
            marker = root / "method.tsv"
            cmd = [sys.executable, str(SCRIPT), "summarize", "--sample", "S1", "--species", "sp", "--ref-name", "ref", "--mt-contig", "chrM", "--regions", str(root / "regions.gz"), "--distribution", str(root / "dist.txt"), "--mt-depth", str(root / "mt.tsv"), "--qc", str(root / "qc.tsv"), "--output", str(root / "mtcn.tsv"), "--marker", str(marker)]
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(marker.exists())


if __name__ == "__main__":
    unittest.main()
