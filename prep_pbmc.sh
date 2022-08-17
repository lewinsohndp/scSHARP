#!/bin/bash
#SBATCH --mem=80g
#SBATCH --time=60
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

python -u pbmc_prep.py
