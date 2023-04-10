#!/bin/bash
#SBATCH --mem=130g
#SBATCH --time=12:00:00
#SBATCH --output=log-%j.out
#SBATCH --error=log-%j.err

python -u cluster_similarity_pbmc.py
