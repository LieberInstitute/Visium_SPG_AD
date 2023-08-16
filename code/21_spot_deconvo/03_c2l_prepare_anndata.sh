#!/bin/bash -l

#$ -cwd
#$ -N "c2l_prepare_anndata"
#$ -o ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#$ -e ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#$ -l mf=40G,h_vmem=40G,h_fsize=50G

#SBATCH -q shared
#SBATCH --mem=40G
#SBATCH --job-name=c2l_prepare_anndata
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log

USE_SLURM=0

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

module load cell2location/0.8a0
python 03_c2l_prepare_anndata.py

echo "**** Job ends ****"
date
