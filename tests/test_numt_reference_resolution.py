from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUMT_NF = ROOT / "numt_detection" / "numt_end2end.nf"
LAUNCHER = ROOT / "launch_pipeline_all_per_batch.sh"
BOUCHET_REF = "/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling"


def test_numt_uses_env_before_nextflow_params():
    text = NUMT_NF.read_text()
    assert "System.getenv('GLOBAL_REF_DIR') ?: params.global_ref_dir" in text
    assert "System.getenv('REF_DIR') ?: params.ref_dir" in text
    assert "System.getenv('NUCLEAR_ONLY_REF_DIR') ?: params.nuclear_only_ref_dir" in text


def test_numt_fails_when_reference_dirs_unconfigured():
    text = NUMT_NF.read_text()
    assert 'error "GLOBAL_REF_DIR/params.global_ref_dir is not configured"' in text
    assert 'error "REF_DIR/params.ref_dir is not configured"' in text
    assert 'error "NUCLEAR_ONLY_REF_DIR/params.nuclear_only_ref_dir is not configured"' in text


def test_launcher_has_no_bouchet_reference_fallback():
    text = LAUNCHER.read_text()
    assert BOUCHET_REF not in text
    assert 'REF_DIR="${REF_DIR:-}"' in text
    assert 'GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-}"' in text
    assert 'NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-}"' in text


def test_auto_numt_config_does_not_write_reference_dirs():
    text = LAUNCHER.read_text()
    start = text.index("write_numt_config()")
    end = text.index("run_numt_nextflow()")
    body = text[start:end]
    assert "SAMPLES_TSV=" in body
    assert "CRAM_ROOT_1=" in body
    assert "DISCOVERY_OUTROOT=" in body
    assert "WHOLE_REF_DIR=" not in body
    assert "NUCLEAR_ONLY_REF_DIR=" not in body
    assert "CHRM_REF_DIR=" not in body


def test_reference_missing_is_prechecked_and_not_retried():
    text = NUMT_NF.read_text()
    assert 'exit 2' in text
    assert 'log_missing_ref "WGS"' in text
    assert 'log_missing_ref "nuclear"' in text
    assert 'CRAM not found' in text
    assert 'CRAI not found' in text
    assert 'mt contig' in text
    retry_block = text[text.index("run_with_retry()"):text.index("# shellcheck disable=SC1090")]
    assert r'"\${status}" -eq 137 || "\${status}" -eq 143' in retry_block
    assert '"${status}" -eq 1' not in retry_block


def test_missing_reference_diagnostics_are_reported():
    text = NUMT_NF.read_text()
    assert "INFO: Candidate filenames checked:" in text
    assert "INFO: Similar files found:" in text
    assert r'find "\${dir}" -maxdepth 1 -iname "*\${name%%_*}*"' in text
