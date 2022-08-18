#!/bin/bash
#SBATCH --mem=8g
#SBATCH --array=1-297%100
#SBATCH --partition=gpu
#SBATCH --time=30
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

real_id=$((${SLURM_ARRAY_TASK_ID} + 1000))
python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/jung/ grid_search_jung_output/ ${real_id}.txt
