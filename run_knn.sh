#!/bin/bash
#SBATCH --mem=150g
#SBATCH --time=8:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

python -u knn_experiment.py
