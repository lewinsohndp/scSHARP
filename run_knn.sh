#!/bin/bash
#SBATCH --mem=250g
#SBATCH --time=24:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

python -u knn_experiment.py
