#!/bin/bash
#SBATCH --mem=16g
#SBATCH --array=1-217%100
#SBATCH --partition=gpu
#SBATCH --time=30
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/jung/ grid_search_jung_output/ ${SLURM_ARRAY_TASK_ID}.txt
