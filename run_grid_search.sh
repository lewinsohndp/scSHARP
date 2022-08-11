#!/bin/bash
#SBATCH --mem=16g
#SBATCH --array=1-540%100
#SBATCH --partition=gpu
#SBATCH --time=30
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

python -u test_config.py /home/groups/ConradLab/daniel/sharp_data/sharp_sims/splat_0.8_de_rq/ grid_search_0.8_output/ ${SLURM_ARRAY_TASK_ID}.txt
