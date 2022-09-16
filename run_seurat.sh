#!/bin/bash
#SBATCH --mem=100g
#SBATCH --time=8:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript prep_umap.R /home/groups/ConradLab/daniel/sharp_data/pbmc_test/counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/
