import shutil
import subprocess
import tempfile
import textwrap
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
VALIDATOR = ROOT / "scripts" / "validate_cram.sh"


class ValidateCramTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.cram = self.root / "sample.cram"
        self.crai = self.root / "sample.cram.crai"
        self.cram.write_bytes(b"C" * 2048)
        self.crai.write_bytes(b"I" * 32)

    def tearDown(self):
        self.tmp.cleanup()

    def mock(self, quickcheck="exit 0"):
        path = self.root / "samtools"
        log = self.root / "samtools.calls"
        path.write_text(textwrap.dedent(f"""\
            #!/bin/bash
            printf '%s\\n' "$*" >> "{log}"
            case "$1" in
              --version) echo 'samtools 1.test'; exit 0;;
              quickcheck) {quickcheck};;
              *) exit 99;;
            esac
        """))
        path.chmod(0o755)
        return path, log

    def run_validator(self, samtools, crai=None, extra=()):
        return subprocess.run(
            [str(VALIDATOR), "--cram", str(self.cram), "--crai", str(crai or self.crai),
             "--samtools", str(samtools), *extra],
            text=True, capture_output=True,
        )

    def assert_status(self, result, returncode, status, reason):
        self.assertEqual(result.returncode, returncode, result.stdout + result.stderr)
        self.assertIn(f"STATUS={status}", result.stdout)
        self.assertIn(f"REASON={reason}", result.stdout)

    def test_complete_marker_is_fastest_path_without_samtools(self):
        marker = self.root / "sample.cram.complete"
        marker.write_text("cram_size=2048\ncrai_size=32\n")
        result = self.run_validator("/samtools/must/not/run", extra=("--marker", str(marker)))
        self.assert_status(result, 0, "COMPLETE", "completion_marker_and_sizes_match")

    def test_legacy_runs_exactly_one_quickcheck(self):
        samtools, log = self.mock()
        result = self.run_validator(samtools)
        self.assert_status(result, 0, "COMPLETE", "files_exist_and_quickcheck_passed")
        calls = log.read_text().splitlines()
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in calls), 1)
        self.assertFalse(any(line.startswith("idxstats") for line in calls))

    def test_ignored_legacy_options_do_not_wait_or_require_reference(self):
        samtools, _ = self.mock()
        started = time.monotonic()
        result = self.run_validator(samtools, extra=(
            "--reference", "/missing/reference.fa", "--stability-retries", "9",
            "--retries", "9", "--delay", "10", "--timeout", "1",
        ))
        self.assert_status(result, 0, "COMPLETE", "files_exist_and_quickcheck_passed")
        self.assertLess(time.monotonic() - started, 2)

    def test_missing_cram_is_incomplete(self):
        self.cram.unlink()
        result = self.run_validator("/samtools/must/not/run")
        self.assert_status(result, 1, "INCOMPLETE", "cram_missing")

    def test_missing_crai_is_incomplete(self):
        result = self.run_validator("/samtools/must/not/run", self.root / "missing.crai")
        self.assert_status(result, 1, "INCOMPLETE", "crai_missing")

    def test_cram_below_minimum_is_incomplete(self):
        self.cram.write_bytes(b"")
        result = self.run_validator("/samtools/must/not/run")
        self.assert_status(result, 1, "INCOMPLETE", "cram_below_minimum_size")

    def test_crai_below_minimum_is_incomplete(self):
        self.crai.write_bytes(b"I")
        result = self.run_validator("/samtools/must/not/run")
        self.assert_status(result, 1, "INCOMPLETE", "crai_below_minimum_size")

    def test_first_quickcheck_failure_is_retried_and_can_complete(self):
        attempts = self.root / "quickcheck.attempts"
        samtools, log = self.mock(textwrap.dedent(f'''\
            if [[ ! -e "{attempts}" ]]; then
              touch "{attempts}"
              echo "I/O error" >&2
              exit 1
            fi
            exit 0
        '''))
        result = self.run_validator(samtools)
        self.assert_status(result, 0, "COMPLETE", "files_exist_and_quickcheck_passed")
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in log.read_text().splitlines()), 2)
        self.assertIn("QUICKCHECK_ATTEMPT=1/2 EXIT=1 OUTPUT=I/O error", result.stderr)
        self.assertIn("QUICKCHECK_ATTEMPT=2/2 EXIT=0 OUTPUT=", result.stderr)

    def test_two_explicit_corruption_failures_are_incomplete(self):
        samtools, log = self.mock('echo "truncated file" >&2; exit 1')
        result = self.run_validator(samtools)
        self.assert_status(result, 1, "INCOMPLETE", "confirmed_corrupt")
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in log.read_text().splitlines()), 2)
        self.assertIn("QUICKCHECK_ATTEMPT=1/2 EXIT=1 OUTPUT=truncated file", result.stderr)
        self.assertIn("QUICKCHECK_ATTEMPT=2/2 EXIT=1 OUTPUT=truncated file", result.stderr)

    def test_two_transient_failures_are_unknown(self):
        samtools, log = self.mock('echo "stale file handle" >&2; exit 1')
        result = self.run_validator(samtools)
        self.assert_status(result, 2, "UNKNOWN", "quickcheck_indeterminate")
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in log.read_text().splitlines()), 2)

    def test_two_empty_failures_are_unknown_and_logged(self):
        samtools, _ = self.mock('exit 3')
        result = self.run_validator(samtools)
        self.assert_status(result, 2, "UNKNOWN", "quickcheck_indeterminate")
        self.assertIn("QUICKCHECK_ATTEMPT=1/2 EXIT=3 OUTPUT=", result.stderr)
        self.assertIn("QUICKCHECK_ATTEMPT=2/2 EXIT=3 OUTPUT=", result.stderr)

    def test_samtools_missing_is_unknown(self):
        result = self.run_validator("/does/not/exist")
        self.assert_status(result, 2, "UNKNOWN", "samtools_unavailable")

    def test_mismatched_marker_falls_back_to_one_quickcheck(self):
        marker = self.root / "bad.complete"
        marker.write_text("cram_size=1\ncrai_size=2\n")
        samtools, log = self.mock()
        result = self.run_validator(samtools, extra=("--marker", str(marker)))
        self.assert_status(result, 0, "COMPLETE", "files_exist_and_quickcheck_passed")
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in log.read_text().splitlines()), 1)

    def test_one_hundred_marked_results_need_no_samtools_and_finish_quickly(self):
        started = time.monotonic()
        for number in range(100):
            cram = self.root / f"sample-{number}.cram"
            crai = self.root / f"sample-{number}.cram.crai"
            marker = self.root / f"sample-{number}.cram.complete"
            cram.write_bytes(b"C" * 1024)
            crai.write_bytes(b"I" * 16)
            marker.write_text("cram_size=1024\ncrai_size=16\n")
            result = subprocess.run([
                str(VALIDATOR), "--cram", str(cram), "--crai", str(crai),
                "--marker", str(marker), "--samtools", "/samtools/must/not/run",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        # This includes 100 separate Python subprocess launches; the validator
        # itself performs no sleep and never starts samtools on the marker path.
        self.assertLess(time.monotonic() - started, 15)

    def test_pipeline_contract_uses_fast_arguments_and_backfills_marker(self):
        source = (ROOT / "preprocessing.nf").read_text()
        command = source[source.index("def command = [validatorScript"):source.index("assert command.every")]
        for obsolete in ("--reference", "--stability-retries", "--retries", "--delay", "--timeout"):
            self.assertNotIn(obsolete, command)
        self.assertIn("StandardCopyOption.ATOMIC_MOVE", source)
        self.assertIn("status == 'UNKNOWN') error", source)
        self.assertIn("status == 'COMPLETE') return false", source)

    def test_bam_to_cram_marker_is_atomic_and_published(self):
        source = (ROOT / "preprocessing.nf").read_text()
        process = source[source.index("process BAM_TO_CRAM"):]
        self.assertLess(process.index('samtools quickcheck -v "${cram_out}"'), process.index('cat > "${meta.id}.cram.complete.tmp"'))
        self.assertIn('mv "${meta.id}.cram.complete.tmp" "${meta.id}.cram.complete"', process)
        self.assertIn('pattern: "*.{cram,crai,complete}"', process)

    @unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not installed")
    def test_nextflow_legacy_result_is_skipped_and_marker_is_backfilled(self):
        sample = "coordinator_test"
        alignment = self.root / "out" / sample / "alignment"
        alignment.mkdir(parents=True)
        (alignment / f"{sample}.cram").write_bytes(b"C" * 2048)
        (alignment / f"{sample}.cram.crai").write_bytes(b"I" * 32)
        samples = self.root / "samples.tsv"
        samples.write_text(f"{sample}\ttest_species\ttest_ref\n")
        samtools, log = self.mock()
        result = subprocess.run([
            shutil.which("nextflow"), "run", str(ROOT / "preprocessing.nf"),
            "--sample_tsv", str(samples), "--outdir", str(self.root / "out"),
            "--samtools_bin", str(samtools),
        ], cwd=self.root, text=True, capture_output=True)
        diagnostics = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, diagnostics)
        self.assertTrue((alignment / f"{sample}.cram.complete").is_file())
        self.assertEqual(sum(line.startswith("quickcheck -v ") for line in log.read_text().splitlines()), 1)
        self.assertIn("status=COMPLETE reason=files_exist_and_quickcheck_passed action=skip", diagnostics)


if __name__ == "__main__":
    unittest.main()
