#!/bin/bash
#SBATCH --mem=130g
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:1

python -u test_model.py
