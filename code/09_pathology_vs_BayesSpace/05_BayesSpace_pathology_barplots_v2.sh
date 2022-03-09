#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=6G,h_vmem=6G,h_fsize=100G
#$ -N BayesSpace_pathology_barplots_v2
#$ -o logs/BayesSpace_pathology_barplots_v2.txt
#$ -e logs/BayesSpace_pathology_barplots_v2.txt
#$ -m e
#$ -hold_jid label_pathology_spots

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
Rscript 05_BayesSpace_pathology_barplots_v2.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
