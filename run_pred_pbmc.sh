#!/bin/bash
#SBATCH --mem=164g
#SBATCH --time=16:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/pbmc_test/counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/pbmc_test/markers.txt /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_meta.csv
