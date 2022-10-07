#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N BayesSpace_enhance_c_merge
#$ -o logs/BayesSpace_enhance_c_merge.$TASK_ID.txt
#$ -e logs/BayesSpace_enhance_c_merge.$TASK_ID.txt
#$ -m e
#$ -t 1-2
#$ -tc 2
#$ -hold_jid BayesSpace_enhance_b_compute_targeted,BayesSpace_enhance_b_compute_wholegenome

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 06_BayesSpace_enhance_c_merge.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
