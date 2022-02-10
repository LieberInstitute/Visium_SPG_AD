#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N BayesSpace_k_search_spe_harmony_wholegenome
#$ -o logs/BayesSpace_k_search_spe_harmony_wholegenome.$TASK_ID.txt
#$ -e logs/BayesSpace_k_search_spe_harmony_wholegenome.$TASK_ID.txt
#$ -m e
#$ -t 4-15
#$ -tc 20
#$ -hold_jid preprocess_and_harmony_spe_postqc

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
Rscript 02_BayesSpace_k_search.R -s spe_harmony_wholegenome

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


