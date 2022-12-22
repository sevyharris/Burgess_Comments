#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=16-00:00:00
#SBATCH --job-name=RMG_2BTP_rerun 
#SBATCH --error=error.rmg_2_BTP.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

source activate rmg_env
python-jl /work/westgroup/nora/Code/RMG-Py/rmg.py input_2_BTP.py

