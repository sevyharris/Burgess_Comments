#!/bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --output=sensitivity.slurm.%x.log
#SBATCH --error=error_sensitivity.slurm.%x.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --mem=20Gb
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1

list_of_input_files=(edited_RMG_with_BROH_777.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_input_files[$index]}"  
/work/westgroup/nora/Code/projects/Burgess_Comments/2_BTP_optimization/models/edited_RMG_with_BROH/
my_path="/work/westgroup/nora/Code/projects/Burgess_Comments/2_BTP_optimization/models/edited_RMG_with_BROH/${folder_name}"

python subtract_sensitivity.py $my_path
