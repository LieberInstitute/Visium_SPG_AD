#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N subset_data_shiny
#$ -o logs/subset_data_shiny.txt
#$ -e logs/subset_data_shiny.txt
#$ -m e
#$ -hold_jid qc_metrics_and_segmentation

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript ../05_deploy_app/subset_data.R
Rscript ../06_deploy_app_targeted/subset_data.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
