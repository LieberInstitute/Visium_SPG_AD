#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N BayesSpace_enhance_a_split
#$ -o logs/BayesSpace_enhance_a_split.$TASK_ID.txt
#$ -e logs/BayesSpace_enhance_a_split.$TASK_ID.txt
#$ -m e
#$ -t 1-2
#$ -tc 2
#$ -hold_jid fasthplus_optimal_k_wholegenome,fasthplus_optimal_k_targeted

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
Rscript 06_BayesSpace_enhance_a_split.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
