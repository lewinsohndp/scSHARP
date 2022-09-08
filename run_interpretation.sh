#!/bin/bash
#SBATCH --mem=100g
#SBATCH --partition=gpu
#SBATCH --time=4:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:1

python -u train_pbmc_model.py
