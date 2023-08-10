#!/bin/bash

#$ -cwd
#$ -N "find_markers"
#$ -o ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#$ -e ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#$ -l mf=40G,h_vmem=40G,h_fsize=50G

#SBATCH -q shared
#SBATCH --mem=40G
#SBATCH --job-name=find_markers
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/02_find_markers.log

USE_SLURM=1

if [[ $USE_SLURM -eq 1 ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 02_find_markers.R

echo "**** Job ends ****"
date
