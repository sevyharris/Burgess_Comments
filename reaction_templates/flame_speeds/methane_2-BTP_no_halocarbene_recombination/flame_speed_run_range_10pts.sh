#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=without_halocarbene_recombination
#SBATCH --error=error.slurm.log
#SBATCH --output=output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_range_10pts.py
