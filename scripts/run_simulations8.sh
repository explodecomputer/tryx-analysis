#!/bin/bash

#SBATCH --job-name=sim7
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=0-10:00:00
#SBATCH --array=1-640
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript run_simulations8.r ${i} 300 1

