#!/bin/bash
#SBATCH --mem=100g
#SBATCH --time=8:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript run_de.R /home/groups/ConradLab/daniel/sharp_data/pbmc_test/
