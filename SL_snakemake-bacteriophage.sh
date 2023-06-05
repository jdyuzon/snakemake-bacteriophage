#!/bin/bash
#SBATCH --account m342
#SBATCH --qos debug
#SBATCH --nodes 2
#SBATCH --constraint cpu
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH -t 00:30:00
#SBATCH -J snakemake-bacteriophage


##conda create -c conda-forge -c bioconda -n snakemake-bacteriophage snakemake
##conda install -n snakemake-bacteriophage -c conda-forge -c bioconda genomad

source /global/homes/j/jdyuzon/.bashrc
conda activate snakemake-bacteriophage

snakemake --unlock
snakemake --cores 32 --rerun-incomplete --forceall -k
