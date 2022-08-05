#!/bin/bash
#SBATCH --mem=64g
#SBATCH --array=1-1000%100
#SBATCH --partition=gpu
#SBATCH --time=1:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:2

python -u test_config.py /home/groups/ConradLab/daniel/sharp_sims/splat_0.5_de_rq/ grid_search_0.5_output/ ${SLURM_ARRAY_TASK_ID}.txt
