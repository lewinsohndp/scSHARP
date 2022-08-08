#!/bin/bash
#SBATCH --mem=64g
#SBATCH --time=120
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_sims/splat_0.6_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.6_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_sims/splat_0.6_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_sims/splat_0.6_de_rq/ref_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.6_de_rq/ref_labels.csv

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_sims/splat_0.7_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.7_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_sims/splat_0.7_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_sims/splat_0.7_de_rq/ref_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.7_de_rq/ref_labels.csv

Rscript tools/r_tools.R /home/groups/ConradLab/daniel/sharp_sims/splat_0.8_de_rq/query_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.8_de_rq/ sctype,scsorter,scina,singler,scpred /home/groups/ConradLab/daniel/sharp_sims/splat_0.8_de_rq/markers.txt /home/groups/ConradLab/daniel/sharp_sims/splat_0.8_de_rq/ref_counts.csv /home/groups/ConradLab/daniel/sharp_sims/splat_0.8_de_rq/ref_labels.csv
