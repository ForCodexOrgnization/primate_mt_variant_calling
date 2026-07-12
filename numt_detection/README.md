# numt_detection
A pipeline for NUMT detection using a pileup strategy.

## One-command end-to-end run

You can now run discovery + downstream best-hit analysis in one step:

```bash
bash run_numt_end2end.sh --config numt_pipeline.config
```

## Config-driven parameters

All input paths and parameters are centralized in a config file.

1. Copy the example config:

```bash
cp numt_pipeline.config.example numt_pipeline.config
```

2. Edit `numt_pipeline.config` with your sample paths and thresholds.

The same config is consumed by:
- `run_numt_discovery.sh`
- `process_numt_candidates_one_ready.sh`
- `run_numt_end2end.sh`

If CRAM/CRAI may exist in two locations, set both primary and fallback paths in config:
- `INPUT_BAM_CRAM` / `INPUT_INDEX`
- `INPUT_BAM_CRAM_ALT` / `INPUT_INDEX_ALT`

When primary files are missing, the pipeline automatically falls back to the alternate paths.

## Run all samples with one Slurm array submission

Sample list must be a 3-column TSV:
1. sample ID
2. real species name
3. reference species name

Then submit once:

```bash
bash submit_numt_end2end_array.sh --config numt_pipeline.config --concurrent 50
```

This will automatically create a Slurm job array and process all rows in `SAMPLES_TSV`.
For each sample, CRAM is searched in `CRAM_ROOT_1` first, then `CRAM_ROOT_2`.


> 说明：`submit_numt_end2end_array.sh` 运行时会按每个样本动态写入 `INPUT_BAM_CRAM/INPUT_INDEX/WGS_REF/NUCLEAR_REF/DISCOVERY_OUTDIR` 到临时 config。
> 所以 array 模式下主 config 不需要填写单样本 CRAM 路径。


> 兼容旧配置：`submit_numt_end2end_array.sh` 仍可读取旧变量 `SAMPLES_FILE`，但建议统一使用 `SAMPLES_TSV`。


## Environment / module load

NUMT runtime environment setup is managed by the repository Nextflow configs,
not by `load_numt_modules.sh`. The deprecated `load_numt_modules.sh` file is
kept as a no-op for backward compatibility only.

For Nextflow-driven NUMT execution, processes labeled `numt_related` use
`params.numt_env_before_script` from `nextflow.config` or `nextflow_mcc.config`
to load SAMtools, BWA, BEDTools, BLAST, Python, and `pysam` via the configured
conda environment. Override `params.numt_env_before_script` in the site-specific
Nextflow config when module names differ.

`run_numt_discovery.sh` validates the runtime environment before running and
reports a clear error if any required tool is missing:

- `samtools`
- `bwa`
- `bedtools`
- `blastn`
- `makeblastdb`
- `python`
- `python -c "import pysam"`
