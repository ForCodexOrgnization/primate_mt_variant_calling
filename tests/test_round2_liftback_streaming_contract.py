from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROUND2 = (ROOT / "primate_pipeline_round2_consensus_NUMT.nf").read_text()
CONFIGS = [ROOT / "nextflow.config", ROOT / "nextflow_mcc.config"]


def liftback_process_source():
    start = ROUND2.index("process LIFTBACK_ROUND2_VCF_TO_ORIGINAL")
    return ROUND2[start:]


def test_liftback_has_dedicated_resources():
    source = liftback_process_source()
    assert "label 'liftback_related'" in source
    assert "label 'generation_related'" not in source.split("script:", 1)[0]

    for config_path in CONFIGS:
        config = config_path.read_text()
        assert "withLabel: liftback_related" in config
        label_config = config.split("withLabel: liftback_related", 1)[1].split("}", 1)[0]
        assert "cpus   = 1" in label_config
        assert "memory = '4.GB'" in label_config
        assert "time   = '2h'" in label_config


def test_ambiguity_metrics_stream_samtools_output():
    source = liftback_process_source()
    metrics = source.split("def metrics(chrom,pos,ref,alt):", 1)[1].split(
        "def append_info", 1
    )[0]

    assert "subprocess.check_output" not in metrics
    assert "out.splitlines()" not in metrics
    assert "subprocess.Popen(command, stdout=subprocess.PIPE, text=True)" in metrics
    assert "for rec in proc.stdout:" in metrics
    assert "rc=proc.wait()" in metrics
    assert "raise subprocess.CalledProcessError(rc, command)" in metrics


def test_ambiguity_info_contract_is_unchanged():
    for name in (
        "alt_support_reads",
        "alt_reads_with_AS",
        "alt_reads_with_XS",
        "alt_reads_with_both_AS_XS",
        "alt_XS_ge_AS_n",
        "alt_XS_ge_AS_pct",
        "alt_XS_eq_AS_n",
        "alt_XS_gt_AS_n",
        "alt_reads_without_XS_ge_AS",
        "AMBIGUITY_METRICS_SOURCE",
    ):
        assert name in liftback_process_source()

    assert "AMBIGUITY_METRICS_SOURCE=not_computed_non_snv" in ROUND2
