#!/bin/bash

#$ -cwd
#$ -N "find_markers"
#$ -o ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#$ -e ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#$ -l caracol,mf=20G,h_vmem=20G,h_fsize=50G

#SBATCH -q shared
#SBATCH --mem=40G
#SBATCH --job-name=find_markers
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/02_find_markers.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/02_find_markers.log

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load conda_R/4.3
Rscript 02_find_markers.R

echo "**** Job ends ****"
date
