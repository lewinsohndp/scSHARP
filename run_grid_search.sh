#!/bin/bash
#SBATCH --mem=32g
#SBATCH --array=1-216
#SBATCH --partition=gpu
#SBATCH --time=60
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:v100:2

# real array is 1-216
#python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/ grid_search_sim_0.7_output/ ${SLURM_ARRAY_TASK_ID}.txt

#python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_de_rq/ grid_search_sim_0.8_output/ ${SLURM_ARRAY_TASK_ID}.txt

#python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/pbmc_test_grid/ grid_search_pbmc_output/ ${SLURM_ARRAY_TASK_ID}.txt

#python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/jung/ grid_search_jung_output_103022/ ${SLURM_ARRAY_TASK_ID}.txt

python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_imbalance_031223/ grid_search_sim_0.7_imbalance_output/ ${SLURM_ARRAY_TASK_ID}.txt

python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.8_imbalance_031223/ grid_search_sim_0.8_imbalance_output/ ${SLURM_ARRAY_TASK_ID}.txt
