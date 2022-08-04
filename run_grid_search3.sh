#!/bin/bash
#SBATCH --mem=64g
#SBATCH --array=1-701%500
#SBATCH --partition=gpu
#SBATCH --time=1:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:1

real_id=$((${SLURM_ARRAY_TASK_ID} + 2000))
SAMPLE=$( sed -n ${real_id}
python -u test_config.py /home/groups/ConradLab/daniel/sharp_sims/splat_0.5_de_rq/ grid_search_0.5_output/ ${real_id}.txt
