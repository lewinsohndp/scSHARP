#!/bin/bash
#SBATCH --mem=16g
#SBATCH --array=75,201
#SBATCH --partition=gpu
#SBATCH --time=60
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

# real array is 1-216
python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/pbmc_test_grid/ grid_search_pbmc_output/ ${SLURM_ARRAY_TASK_ID}.txt
