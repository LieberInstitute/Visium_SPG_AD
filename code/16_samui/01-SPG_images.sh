#!/bin/bash

#SBATCH -p shared
#SBATCH --array=1%10
#SBATCH --mem=20G
#SBATCH --job-name=SPG_images
#SBATCH -o ../../processed-data/16_samui/logs/01-SPG_images_%a.log
#SBATCH -e ../../processed-data/16_samui/logs/01-SPG_images_%a.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load samui/1.0.0-next.24
python 01-SPG_images.py

echo "**** Job ends ****"
date
