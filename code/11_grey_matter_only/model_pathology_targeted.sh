#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -N model_pathology_targeted
#$ -o logs/model_pathology_targeted.txt
#$ -e logs/model_pathology_targeted.txt
#$ -m e
#$ -hold_jid explore_expr_variability_targeted

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
Rscript 03_model_pathology.R -s targeted

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


