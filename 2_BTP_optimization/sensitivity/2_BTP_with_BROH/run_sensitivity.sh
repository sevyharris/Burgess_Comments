#!/bin/bash
#SBATCH --job-name=sens_NIST
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1



python sensitivity.py 
