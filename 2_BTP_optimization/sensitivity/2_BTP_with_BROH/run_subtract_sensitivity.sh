#!/bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --output=sensitivity.slurm.%x.log
#SBATCH --error=error_sensitivity.slurm.%x.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --mem=20Gb
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-2

list_of_input_files=(RMG_NK/chem_annotated.cti RMG_with_BROH/chemkin/copies/copy_chem0136.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_input_files[$index]}"  

my_path="/work/westgroup/nora/Code/projects/Burgess_Comments/2_BTP_optimization/models/${folder_name}"


python subtract_sensitivity.py $my_path
