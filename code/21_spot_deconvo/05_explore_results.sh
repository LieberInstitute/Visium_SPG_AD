#!/bin/bash

#$ -cwd
#$ -N "explore_results"
#$ -o ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#$ -e ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#$ -l mf=5G,h_vmem=5G

#SBATCH -q shared
#SBATCH --mem=5G
#SBATCH --job-name=explore_results
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/05_explore_results.log

#   Auto-detect if SGE or SLURM is being used
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
Rscript 05_explore_results.R

echo "**** Job ends ****"
date
