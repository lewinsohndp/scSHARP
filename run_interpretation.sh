#!/bin/bash
#SBATCH --mem=160g
#SBATCH --partition=gpu
#SBATCH --time=16:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:v100:1

python -u train_pbmc_model.py
