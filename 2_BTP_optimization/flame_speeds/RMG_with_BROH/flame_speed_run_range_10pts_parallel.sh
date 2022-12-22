#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=6-00:00:00
#SBATCH --job-name=edited_RMG_ss
#SBATCH --error=error.slurm.%x.log
#SBATCH --output=output.slurm.%x.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=20Gb 
#SBATCH --array=1-20
#SBATCH --partition=west

list_of_input_files=(copy_chem0130.cti copy_chem0131.cti copy_chem0132.cti copy_chem0133.cti copy_chem0134.cti copy_chem0135.cti copy_chem0136.cti copy_chem0137.cti copy_chem0138.cti copy_chem0139.cti copy_chem0140.cti copy_chem0141.cti copy_chem0142.cti copy_chem0143.cti copy_chem0144.cti copy_chem0145.cti copy_chem0146.cti copy_chem0147.cti copy_chem0148.cti copy_chem0149.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_input_files[$index]}"  

my_path="/work/westgroup/nora/Code/projects/Burgess_Comments/2_BTP_optimization/models/RMG_with_BROH/chemkin/copies/${folder_name}"


source activate cantera_env
python flame_speed_calc_range_10pts_parallel.py $my_path
