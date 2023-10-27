#!/bin/bash
#SBATCH --job-name=snakemake_SEAAD     # Job name
#SBATCH --partition=1gpu               # Add the appropriate partition name

snakemake --unlock
snakemake -j 12 -R     # Run snakemake. Adjust as necessary.
