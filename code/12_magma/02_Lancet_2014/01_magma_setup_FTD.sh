#!/bin/bash
#$ -cwd
#$ -N magma_
#$ -o ./logs/magma-gsa_steps1-3_ADHD_MNT06Jun2022.o
#$ -e ./logs/magma-gsa_steps1-3_ADHD_MNT06Jun2022.e
#$ -l bluejay,mem_free=16G,h_vmem=20G

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
Rscript 01_magma_setup_FTD.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
