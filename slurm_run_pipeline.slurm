#!/bin/bash
#SBATCH -J snakemake
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -o %N.%J.out # STDOUT
#SBATCH -e %N.%j.err # STDERR
#SBATCH --time=168:00:00

eval "$(conda shell.bash hook)"
conda activate snakemake

# Pass the same arguments (e.g. --configfile) so the right output_dir is unlocked
snakemake -j 1 --unlock --cleanup-shadow "$@"
snakemake --profile config/snakemake_profile/slurm "$@"