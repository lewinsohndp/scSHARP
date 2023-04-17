#!/bin/bash
#SBATCH --mem=200g
#SBATCH --time=8:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/pbmc_test/counts_copy_mod_3.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/pbmc_test/markers.txt /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_meta.csv

