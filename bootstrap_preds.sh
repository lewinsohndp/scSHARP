#!/bin/bash
#SBATCH --mem=16g
#SBATCH --array=1-10
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq_boot/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq/ref_counts.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq/ref_labels.csv

mv /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq_boot/preds.csv /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.7_de_rq_boot/preds${SLURM_ARRAY_TASK_ID}.csv
