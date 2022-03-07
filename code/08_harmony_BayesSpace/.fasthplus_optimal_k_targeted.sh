#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N fasthplus_optimal_k_targeted
#$ -o logs/fasthplus_optimal_k_targeted.$TASK_ID.txt
#$ -e logs/fasthplus_optimal_k_targeted.$TASK_ID.txt
#$ -m e
#$ -t 2-15
#$ -tc 14
#$ -hold_jid BayesSpace_k_search_wholegenome,BayesSpace_k_search_targeted

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
Rscript 05_fasthplus_optimal_k -s targeted

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


