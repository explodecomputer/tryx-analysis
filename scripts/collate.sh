#!/bin/bash

#SBATCH --job-name=sim6
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=0-10:00:00
#SBATCH --output=job_reports/slurm-collate.out
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript collate_simulations6.r
Rscript collate_simulations7.r


