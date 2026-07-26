import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


VALIDATOR = Path(__file__).parents[1] / "scripts" / "validate_cram.sh"


class ValidateCramTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.cram = root / "sample.cram"
        self.crai = root / "sample.cram.crai"
        self.ref = root / "ref.fa"
        self.cram.write_bytes(b"C" * 2048)
        self.crai.write_bytes(b"I" * 32)
        self.ref.write_text(">chrM\nA\n")
        Path(str(self.ref) + ".fai").write_text("chrM\t1\t6\t1\t2\n")

    def tearDown(self):
        self.tmp.cleanup()

    def mock(self, quickcheck="exit 0", idxstats="exit 0"):
        path = Path(self.tmp.name) / "samtools"
        path.write_text(textwrap.dedent(f"""\
            #!/bin/bash
            case "$1" in
              --version) echo 'samtools 1.test'; exit 0;;
              quickcheck) {quickcheck};;
              idxstats) {idxstats};;
            esac
        """))
        path.chmod(0o755)
        return path

    def run_validator(self, samtools, crai=None, extra=()):
        return subprocess.run([str(VALIDATOR), "--cram", str(self.cram), "--crai", str(crai or self.crai),
            "--reference", str(self.ref), "--samtools", str(samtools), "--retries", "3", "--delay", "0",
            *extra], text=True, capture_output=True)

    def test_complete_with_marker(self):
        marker = Path(self.tmp.name) / "sample.cram.complete"
        marker.write_text("cram_size=2048\ncrai_size=32\n")
        result = self.run_validator(self.mock(), extra=("--marker", str(marker)))
        self.assertEqual(result.returncode, 0); self.assertIn("STATUS=COMPLETE", result.stdout)

    def test_legacy_without_marker(self):
        self.assertEqual(self.run_validator(self.mock()).returncode, 0)

    def test_samtools_missing_is_unknown(self):
        result = self.run_validator("/does/not/exist")
        self.assertEqual(result.returncode, 2); self.assertIn("samtools_unavailable", result.stdout)

    def test_transient_then_success(self):
        count = Path(self.tmp.name) / "count"
        cmd = f'n=$(cat "{count}" 2>/dev/null || echo 0); n=$((n+1)); echo $n > "{count}"; if ((n==1)); then echo "Input/output error" >&2; exit 1; fi; exit 0'
        result = self.run_validator(self.mock(cmd))
        self.assertEqual(result.returncode, 0); self.assertIn("QUICKCHECK_ATTEMPTS=2", result.stdout)

    def test_confirmed_truncation_is_incomplete(self):
        result = self.run_validator(self.mock('echo "truncated file" >&2; exit 1'))
        self.assertEqual(result.returncode, 1); self.assertIn("confirmed_corrupt", result.stdout)

    def test_missing_crai_is_incomplete(self):
        result = self.run_validator(self.mock(), Path(self.tmp.name) / "missing.crai")
        self.assertEqual(result.returncode, 1); self.assertIn("crai_missing", result.stdout)

    def test_alternate_crai_is_valid(self):
        alternate = Path(self.tmp.name) / "sample.crai"; alternate.write_bytes(b"I" * 32)
        self.assertEqual(self.run_validator(self.mock(), alternate).returncode, 0)

    def test_timeout_is_unknown(self):
        result = self.run_validator(self.mock("sleep 2; exit 0"), extra=("--timeout", "1"))
        self.assertEqual(result.returncode, 2); self.assertIn("STATUS=UNKNOWN", result.stdout)

    def test_marker_mismatch_is_unknown(self):
        marker = Path(self.tmp.name) / "bad.complete"; marker.write_text("cram_size=1\ncrai_size=2\n")
        result = self.run_validator(self.mock(), extra=("--marker", str(marker)))
        self.assertEqual(result.returncode, 2); self.assertIn("marker_size_mismatch", result.stdout)

    def test_force_reprocess_guard_precedes_validation(self):
        source = (Path(__file__).parents[1] / "preprocessing.nf").read_text()
        force = source.index("if (paramAsBoolean(params.force_reprocess_existing_cram))")
        validation = source.index("def validation = validateExistingCram", force)
        self.assertLess(force, validation)
        self.assertIn("reason=user_forced", source[force:validation])

    def test_growing_cram_is_unknown(self):
        grower = subprocess.Popen(["bash", "-c", f'for i in {{1..15}}; do printf X >> "{self.cram}"; sleep .2; done'])
        try:
            result = self.run_validator(self.mock(), extra=("--delay", "1", "--stability-retries", "2"))
        finally:
            grower.wait()
        self.assertEqual(result.returncode, 2); self.assertIn("files_not_stable", result.stdout)


if __name__ == "__main__":
    unittest.main()
