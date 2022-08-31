#!/bin/bash
#SBATCH --mem=16g
#SBATCH --array=1-216
#SBATCH --partition=gpu
#SBATCH --time=60
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

# real array is 1-216
python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_0.7_de_rq/ grid_search_sim_output/ ${SLURM_ARRAY_TASK_ID}.txt
