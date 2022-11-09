#!/bin/bash
#SBATCH --mem=120g
#SBATCH --partition=gpu
#SBATCH --time=2:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err
#SBATCH --gres=gpu:1

#--gres=gpu:v100:1
python -u train_pbmc_model.py
