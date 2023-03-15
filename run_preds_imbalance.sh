#!/bin/bash
#SBATCH --mem=64g
#SBATCH --time=4:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err


Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/query_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/markers.txt /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/ref_counts_filtered.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/ref_labels_filtered.csv

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/query_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/markers.txt /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/ref_counts_filtered.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/ref_labels_filtered.csv


