#!/bin/bash
#SBATCH --job-name=encDx
#SBATCH --output=/mnt/arcus/lab/users/kafadare/out_messages/%x_%A_%a_output.txt
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

CHUNKS=$1

module load R/4.2.3

# Run the R script with the specified model
Rscript /mnt/arcus/lab/users/kafadare/encDx_inparts.R $SLURM_ARRAY_TASK_ID $CHUNKS