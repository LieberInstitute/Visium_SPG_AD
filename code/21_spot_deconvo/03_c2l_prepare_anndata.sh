#!/bin/bash -l

#$ -cwd
#$ -N "c2l_prepare_anndata"
#$ -o ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#$ -e ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#$ -l caracol,mf=40G,h_vmem=40G,h_fsize=50G

#SBATCH -q shared
#SBATCH --mem=40G
#SBATCH --job-name=c2l_prepare_anndata
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/03_c2l_prepare_anndata.log

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

module load cell2location/0.1.3
python 03_c2l_prepare_anndata.py

echo "**** Job ends ****"
date
