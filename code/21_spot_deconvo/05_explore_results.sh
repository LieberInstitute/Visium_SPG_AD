#!/bin/bash -l

#$ -cwd
#$ -N "explore_results"
#$ -o ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#$ -e ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#$ -l mf=5G,h_vmem=5G
#$ -t 1-3
#$ -tc 3

#SBATCH -q shared
#SBATCH --mem=5G
#SBATCH --job-name=explore_results
#SBATCH -o ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#SBATCH -e ../../processed-data/21_spot_deconvo/logs/05_explore_results.log
#SBATCH --array=1-3%3

#   Auto-detect if SGE or SLURM is being used
if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
    task_id=$SLURM_ARRAY_TASK_ID
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
    task_id=$SGE_TASK_ID
fi

all_groups=(white gray both)
subset=${all_groups[$(($task_id - 1))]}

log_path="../../processed-data/21_spot_deconvo/logs/05_explore_results_${subset}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"
echo "Task id: ${task_id}"

module load conda_R/4.3
Rscript 05_explore_results.R -s $subset

echo "**** Job ends ****"
date
} > $log_path 2>&1
