# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Snakemake pipeline that trains ChromBPNet models on ATAC/DNase BAMs, then runs prediction, contribution scores, TF-MoDISco motif discovery and TOBIAS BINDetect footprinting on the contribution tracks. Runs on a SLURM cluster (IRB Barcelona); GPU rules target the `irb_gpu_3090` partition.

## Running

All commands run from the repo root (paths like `config/...`, `data/...` are relative to it).

- Full run (SLURM head job; extra args pass through to snakemake): `sbatch slurm_run_pipeline.slurm`
  - It activates the `snakemake` conda env, unlocks and cleans shadow dirs, then runs `snakemake --profile config/snakemake_profile/slurm "$@"`.
- Dry run / DAG check (safe in the login shell): `snakemake -n --profile config/snakemake_profile/slurm` (add the same `--configfile` you would pass to the head job).
- Single target: `snakemake -n --profile config/snakemake_profile/slurm <output_path>`, e.g. `pred_bw/<sample>_fold_0_chrombpnet_nobias.bw` (relative to `output_dir`). Swap `-n` for a real run only through the sbatch head job.
- Small end-to-end test (2 downsampled samples, 2 folds, 2 epochs, generated splits): `sbatch slurm_run_pipeline.slurm --configfile config/test/config.yaml`. Output dir `/scratch/qscd/dfernandezp/chrombpnetTest`; its downsampled BAMs live in `input_bams/` there (see `config/test/samples.tsv`).
- GPU sanity check: request target `logs/gpu_test/tensorflow_gpu_check_successful.flag` (rule `gpu_test_rule`, script `workflow/scripts/test_gpu_chrombpnet.py`).

There are no unit tests; validate changes with `snakemake -n -p` (and `--lint`). The fold-averaging scripts in `workflow/scripts/average_*.py` are standalone CLIs and can be tested on tiny files inside the chrombpnet image (`singularity exec ~/chrombpnet_latest.sif python ...`).

The profile (`config/snakemake_profile/slurm/config.yaml`) uses the SLURM executor, `use-singularity` with `--bind /home --nv` (GPU passthrough), conda for TOBIAS, `keep-going`, and default account `qscd` / partition `irb_cpu_iclk`.

## Configuration

- `config/config.yaml`: genome build, output dir (`/scratch/qscd/dfernandezp/chrombpnet`), folds to run (`folds: ["0"]` by default), MACS3/ChromBPNet/MoDISco params, container images. Note the keys are spelled `chromnpnet_bias` / `chromnpnet_model` (typo is load-bearing; rules reference these exact names).
- `config/samples.tsv`: `SampleName`, `BamPath` (multiple BAMs separated by `;`, merged per sample), `UseForBiasModel` (`True` samples' BAMs are all pooled into a single bias model).
- `config/resources.yaml`: per-rule `mem_mb`/`runtime` (minutes)/`cpu`/`retries`, plus `gres`/`slurm_partition` for GPU rules. Every rule reads its resources via `RESOURCES["<rule_name>"]`, so adding a rule means adding an entry here.
- Fold splits: `splits_dir` (default `config/chr_splits`, pre-made mm10 `fold_N.json`) is used as is. If `splits_dir: ""`, rules `get_canonical_chrom_sizes` + `make_chr_splits` generate them with `chrombpnet prep splits` from `chr_splits.folds` (test/valid chroms per fold) into `<output_dir>/chr_splits/`. Rules refer to the fold file through the `FOLD_JSON` global.
- `config/Blacklist/`: vendored ENCODE Blacklist repo; the pipeline only uses `lists/<build>-blacklist.v2.bed.gz` unless `custom_blacklist` is set.

## Architecture

`output_dir` holds everything that is not pipeline code/config: the Snakefile sets `workdir: OUTPUT_DIR`, so `.snakemake/` (singularity images, conda envs, metadata, slurm logs, shadow dirs), downloaded reference data (`data/`), the singularity pull cache (`.singularity_cache/`), the conda package cache (`.conda_pkgs/`), tmp files and results all go there. Point it to `/scratch` for real runs (home has a small quota). Consequences: relative paths in the config (sample sheet, its BAMs, resources, blacklist, `splits_dir`, `output_dir`) are relative to the repo root and must be wrapped in `repo_path()` in the Snakefile *before* the `workdir:` line (SLURM rule jobs re-parse the Snakefile from inside `output_dir`, so the current directory can't be used); rule inputs/outputs are written as plain relative paths (no `OUTPUT_DIR` prefix, so no f-strings unless a value like `TOP_N_PEAKS` is interpolated) and resolve inside `output_dir`. Singularity binds `/home` and `/scratch/qscd` (never the whole `/scratch`: the chrombpnet image has its code and data under `/scratch/chrombpnet`, which would be hidden), so `output_dir` must be under one of them.

`workflow/Snakefile` loads config, parses the sample sheet into `SAMPLES` and `ALL_BIAS_BAMS`, defines shared globals (`OUTPUT_DIR`, `GENOME_FASTA`, `TOP_N_PEAKS`, `FOLDS`, ...) used by all rule files, sets `TMPDIR`/`SINGULARITY_TMPDIR` to `<output_dir>/<host_tmp_dir>`, then includes rules in order (later files reference `rules.<name>.output` from earlier ones):

1. `prepare_data.smk`: downloads genome FASTA, chrom sizes, JASPAR motifs into `<output_dir>/data/`; unzips the blacklist; merges/sorts/indexes per-sample BAMs.
2. `generate_bias_model.smk`: pools all bias BAMs → MACS3 → blacklist filter (blacklist slopped by 1057 bp = half ChromBPNet input) → top-N peaks (`scripts/filter_top_peaks.py`, by column 7) → `chrombpnet prep nonpeaks` → `chrombpnet bias pipeline` per fold.
3. `preprocessing.smk`: the same peak/nonpeak steps per sample.
4. `training.smk`: `chrombpnet pipeline` per sample × fold, using the fold's shared bias model. Produces `*_chrombpnet_nobias.h5`.
5. `prediction.smk`: per sample × fold `pred_bw`, `contribs_bw` (profile/counts score h5 + bigwigs) and `chrombpnet footprints` (motifs from `footprints.motifs_to_pwm` inserted in the fold's test-chromosome nonpeaks); per sample `modisco motifs` + `modisco report` on the fold-averaged profile scores.
6. `fold_average.smk`: per sample, averages folds (`{sample}_mean*` outputs) with `scripts/average_bigwigs.py` (prediction and contribution bigwigs) and `scripts/average_contribs_h5.py` (contribution h5; requires identical `interpreted_regions.bed` across folds). With one fold the mean equals that fold.
7. `tobias.smk`: merges all samples' (fold-averaged) interpreted regions into consensus peaks, then runs TOBIAS BINDetect twice (profile and counts contribution bigwigs as "footprint" signals, one condition per sample). Uses the conda env `workflow/envs/tobias.yaml` (bioconda `tobias` package; installing it with pip no longer builds), not a container, and forces `TMPDIR=/tmp`.

`wildcard_constraints` in the Snakefile restrict `sample` to the sample names and `fold` to digits, which keeps `{sample}_fold_{fold}*` and `{sample}_mean*` patterns unambiguous.

Conventions to preserve:
- ChromBPNet refuses to write into an existing output dir, but Snakemake pre-creates output parents. Training rules therefore write to a `tmp_*` dir, `mv` it to the final location, and touch a `complete.flag`. Their `evaluation/*overall_report.html` QC reports are declared outputs (used to tune `bias_thresh_factor`).
- `chrombpnet prep nonpeaks` rules use `shadow: "minimal"` (hence `--cleanup-shadow` in the head script).
- Most rules run inside `chrombpnet_container` (which also provides samtools/bedtools/modisco); MACS3 and bedtools-only steps use their own biocontainers. The MoDISco MEME db path points inside the ChromBPNet image.
- `/scratch` purges files with mtimes older than ~60 days, and conda keeps the packages' original mtimes, so an `onstart:` handler in the Snakefile touches `.snakemake/conda` and `.conda_pkgs` on every run. If an env breaks anyway (e.g. `ModuleNotFoundError` inside it after a long pause), delete those two dirs in `output_dir`.
- The Snakefile uses nested same-quote f-strings, which requires Python ≥ 3.12 in the snakemake env.
