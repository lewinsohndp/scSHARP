#!/bin/bash
#SBATCH --mem=250g
#SBATCH --time=24:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

#Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_sims_test/splat_0.6_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_sims_test/splat_0.6_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_sims_test/splat_0.6_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_sims_test/splat_0.6_de_rq/ref_counts.csv /home/groups/ConradLab/daniel/sharp_sims_test/splat_0.6_de_rq/ref_labels.csv

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/ref_counts_filtered.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/ref_labels_filtered.csv

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/ref_counts_filtered.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/ref_labels_filtered.csv

#Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/jung/counts.csv /home/groups/ConradLab/daniel/sharp_data/jung/ sctype,scsorter,scina /home/groups/ConradLab/daniel/sharp_data/jung/markers.txt

#Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/pbmc_test/counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/pbmc_test/markers.txt /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_counts.csv /home/groups/ConradLab/daniel/sharp_data/pbmc_test/pbmc3k_ref_meta.csv


