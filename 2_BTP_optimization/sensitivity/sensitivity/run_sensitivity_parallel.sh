#!/bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --output=fc_all.slurm.%x.log
#SBATCH --error=error_fc_all.slurm.%x.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-2

list_of_input_files=(2-BTP_kinetics_originally_no_M.cti 2-BTP_kinetics_with_M.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_input_files[$index]}"  

my_path="/work/westgroup/nora/Code/projects/Burgess_Comments/2_BTP_optimization/models/NIST/${folder_name}"


python sensitivity_parallel.py $my_path
