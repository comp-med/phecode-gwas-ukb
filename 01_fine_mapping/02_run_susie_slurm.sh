#!/bin/sh

## script to runs econd step of REGENIE

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=72:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task 120

#! no other jobs can be run on the same node
#SBATCH --exclusive

#! run as job array
#SBATCH --array=1%4

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
expo="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/phecode.fine.mapping.regions.large.txt)"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/phecode.fine.mapping.regions.large.txt)"
low="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' input/phecode.fine.mapping.regions.large.txt)"
upp="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' input/phecode.fine.mapping.regions.large.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Exposure: ${expo} | Chromosome : ${chr} | Start : ${low} | End : ${upp}"

echo "${chr}:${low}-${upp}"

## set up R environment
source activate Renv

## run the R script
scripts/03_susie_fine_mapping.R ${expo} ${chr} ${low} ${upp} 

## deactivate conda env
conda deactivate

## do some cleaning
rm tmpdir/tmp.${expo}.${chr}.${low}.${upp}.dosage