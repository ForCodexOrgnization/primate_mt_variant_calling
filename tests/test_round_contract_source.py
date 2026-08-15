from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
round1 = (ROOT / "primate_pipeline_numt_decoy_round1.nf").read_text()
round2 = (ROOT / "primate_pipeline_round2_consensus_NUMT.nf").read_text()
config = (ROOT / "nextflow.config").read_text()
launcher = (ROOT / "launch_pipeline_streaming_per_sample.sh").read_text()

# The current producer, shared completion checker configuration, consumer, and
# final-coverage lookup all agree on the modern tree.
assert 'round_1/mtdna_variant_calling' in round1
assert 'round1_vcf_subdir = "round_1/mtdna_variant_calling"' in config
assert 'round_1/mtdna_variant_calling' in launcher
assert 'VCF_ROOT="\\${SAMPLE_DIR}/${params.round1_vcf_subdir}"' in round2

# mtDNA WDL calls and decoy-coordinate NUMT consensus calls remain explicitly
# separate, with the same filtered PASS NUMT product used by the Round1 builder.
assert "ROUND1_MTDNA_CONSENSUS_VCF" in round2
assert "ROUND1_NUMT_CONSENSUS_VCF" in round2
assert ".numt_decoy.pass.split.vcf.gz" in round2
assert "Precomputed consensus NUMT FASTA is already selected" in round2

# Compatibility is current-first, then legacy; it is never a modern prerequisite.
current = round2.index('VCF_ROOT="\\${SAMPLE_DIR}/${params.round1_vcf_subdir}"')
legacy = round2.index('LEGACY_VCF_ROOT="\\${SAMPLE_DIR}/round_1_variant_calling_decoy"')
find_current = round2.index('find "\\${VCF_ROOT}"', legacy)
fallback = round2.index("searching legacy fallback", find_current)
assert current < legacy < find_current < fallback
