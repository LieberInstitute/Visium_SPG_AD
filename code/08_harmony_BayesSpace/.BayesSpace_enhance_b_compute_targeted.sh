#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N BayesSpace_enhance_b_compute_targeted
#$ -o logs/BayesSpace_enhance_b_compute_targeted.$TASK_ID.txt
#$ -e logs/BayesSpace_enhance_b_compute_targeted.$TASK_ID.txt
#$ -m e
#$ -t 1-10
#$ -tc 10
#$ -hold_jid BayesSpace_enhance_a_split

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
Rscript 06_BayesSpace_enhance_b_compute.R -s targeted

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


