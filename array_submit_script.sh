#!/bin/bash

SUBMIT_SCRIPT=$1
ARRAY_SIZE=$2


sbatch --array=1-$ARRAY_SIZE "$SUBMIT_SCRIPT" "$ARRAY_SIZE"
echo "Job Array ID: $job_array_id"

sleep 1

#Waiting for all jobs to finish to concatanate output files
while squeue --job "$job_array_id" | grep -q "R\| PD"; do
  echo "Waiting for job array $job_array_id to finish."
  sleep 400
done

echo "All jobs in array have completed."