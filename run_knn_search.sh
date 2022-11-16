#!/bin/bash
#SBATCH --mem=150g
#SBATCH --time=24:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

python -u knn_search.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/

python -u knn_search.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/

#python -u knn_search.py /home/groups/ConradLab/daniel/sharp_data/pbmc_test_grid/

#python -u knn_search.py /home/groups/ConradLab/daniel/sharp_data/jung/
