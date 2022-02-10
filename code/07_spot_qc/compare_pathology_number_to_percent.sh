#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N compare_pathology_number_to_percent
#$ -o logs/compare_pathology_number_to_percent.txt
#$ -e logs/compare_pathology_number_to_percent.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.1.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript compare_pathology_number_to_percent.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
