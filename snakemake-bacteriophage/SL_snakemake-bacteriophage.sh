#!/bin/bash
#SBATCH --account m342
#SBATCH --qos debug
#SBATCH --nodes 2
#SBATCH --constraint cpu
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH -t 00:10:00
#SBATCH -J snakemake-bacteriophage


source /global/homes/j/jdyuzon/.bashrc
conda activate snakemake-bacteriophage

snakemake --unlock
snakemake --cores 128 --rerun-incomplete --forceall
