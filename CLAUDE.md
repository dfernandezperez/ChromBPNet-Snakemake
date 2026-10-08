# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Snakemake pipeline that trains ChromBPNet models on ATAC/DNase BAMs, then runs prediction, contribution scores, TF-MoDISco motif discovery and TOBIAS BINDetect footprinting on the contribution tracks. Runs on a SLURM cluster (IRB Barcelona); GPU rules target the `irb_gpu_3090` partition.

## Running

All commands run from the repo root (paths like `config/...`, `data/...` are relative to it).

- Full run (SLURM head job; extra args pass through to snakemake): `sbatch slurm_run_pipeline.slurm`
  - It activates the `snakemake` conda env, unlocks and cleans shadow dirs, then runs `snakemake --profile config/snakemake_profile/slurm "$@"`.
- Dry run / DAG check (safe in the login shell): `snakemake -n --profile config/snakemake_profile/slurm`
- Single target: `snakemake -n --profile config/snakemake_profile/slurm <output_path>`, e.g. `results/pred_bw/<sample>_fold_0_chrombpnet_nobias.bw`. Swap `-n` for a real run only through the sbatch head job.
- GPU sanity check: request target `results/logs/gpu_test/tensorflow_gpu_check_successful.flag` (rule `gpu_test_rule`, script `workflow/scripts/test_gpu_chrombpnet.py`).

There are no unit tests or linters; validate changes with `snakemake -n` (and `--lint` if needed).

The profile (`config/snakemake_profile/slurm/config.yaml`) uses the SLURM executor, `use-singularity` with `--bind /home --nv` (GPU passthrough), conda for TOBIAS, `keep-going`, and default account `qscd` / partition `irb_cpu_iclk`.

## Configuration

- `config/config.yaml`: genome build, output dir (`results`), folds to run (`folds: ["0"]` by default), MACS3/ChromBPNet/MoDISco params, container images. Note the keys are spelled `chromnpnet_bias` / `chromnpnet_model` (typo is load-bearing; rules reference these exact names).
- `config/samples.tsv`: `SampleName`, `BamPath` (multiple BAMs separated by `;`, merged per sample), `UseForBiasModel` (`True` samples' BAMs are all pooled into a single bias model).
- `config/resources.yaml`: per-rule `mem_mb`/`runtime` (minutes)/`cpu`/`retries`, plus `gres`/`slurm_partition` for GPU rules. Every rule reads its resources via `RESOURCES["<rule_name>"]`, so adding a rule means adding an entry here.
- `config/chr_splits/fold_N.json`: train/valid/test chromosome splits (mm10 chromosomes).
- `config/Blacklist/`: vendored ENCODE Blacklist repo; the pipeline only uses `lists/<build>-blacklist.v2.bed.gz` unless `custom_blacklist` is set.

## Architecture

`workflow/Snakefile` loads config, parses the sample sheet into `SAMPLES` and `ALL_BIAS_BAMS`, defines shared globals (`OUTPUT_DIR`, `GENOME_FASTA`, `TOP_N_PEAKS`, `FOLDS`, ...) used by all rule files, sets `TMPDIR`/`SINGULARITY_TMPDIR` to `<output_dir>/<host_tmp_dir>`, then includes rules in order (later files reference `rules.<name>.output` from earlier ones):

1. `prepare_data.smk`: downloads genome FASTA, chrom sizes, JASPAR motifs into `data/` (outside `results/`, gitignored); unzips the blacklist; merges/sorts/indexes per-sample BAMs.
2. `generate_bias_model.smk`: pools all bias BAMs → MACS3 → blacklist filter (blacklist slopped by 1057 bp = half ChromBPNet input) → top-N peaks (`scripts/filter_top_peaks.py`, by column 7) → `chrombpnet prep nonpeaks` → `chrombpnet bias pipeline` per fold.
3. `preprocessing.smk`: the same peak/nonpeak steps per sample.
4. `training.smk`: `chrombpnet pipeline` per sample × fold, using the fold's shared bias model. Produces `*_chrombpnet_nobias.h5`.
5. `prediction.smk`: `pred_bw`, `contribs_bw` (profile/counts score h5 + bigwigs), `modisco motifs` + `modisco report`.
6. `tobias.smk`: merges all samples' interpreted regions into consensus peaks, then runs TOBIAS BINDetect twice (profile and counts contribution bigwigs as "footprint" signals, one condition per sample). Uses the conda env `workflow/envs/tobias.yaml`, not a container, and forces `TMPDIR=/tmp`.

Conventions to preserve:
- ChromBPNet refuses to write into an existing output dir, but Snakemake pre-creates output parents. Training rules therefore write to a `tmp_*` dir, `mv` it to the final location, and touch a `complete.flag`.
- `chrombpnet prep nonpeaks` rules use `shadow: "minimal"` (hence `--cleanup-shadow` in the head script).
- Most rules run inside `chrombpnet_container` (which also provides samtools/bedtools/modisco); MACS3 and bedtools-only steps use their own biocontainers. The MoDISco MEME db path points inside the ChromBPNet image.
- The Snakefile uses nested same-quote f-strings, which requires Python ≥ 3.12 in the snakemake env.
