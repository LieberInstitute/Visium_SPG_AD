#!/bin/bash

#$ -cwd
#$ -N "read_data"
#$ -o ../../processed-data/21_spot_deconvo/logs/01_read_data.log
#$ -e ../../processed-data/21_spot_deconvo/logs/01_read_data.log
#$ -pe local 4
#$ -l mf=10G,h_vmem=10G,h_fsize=50G

#SBATCH -q shared
#SBATCH --mem=40G
#SBATCH -c 4
#SBATCH --job-name=read_data
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/01_read_data.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/01_read_data.log

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
Rscript 01_read_data.R

echo "**** Job ends ****"
date
