#!/bin/bash
#$ -cwd
#$ -l bleujay,mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N build_basic_spe
#$ -o logs/build_basic_spe.txt
#$ -e logs/build_basic_spe.txt
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
Rscript build_basic_spe.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
