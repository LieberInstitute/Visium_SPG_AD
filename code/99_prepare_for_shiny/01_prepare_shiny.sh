#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=7G,h_vmem=7G,h_fsize=100G
#$ -N prepare_shiny
#$ -o logs/prepare_shiny.$TASK_ID.txt
#$ -e logs/prepare_shiny.$TASK_ID.txt
#$ -m e
#$ -hold_jid prepare_to_share
#$ -t 1-2
#$ -tc 2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_prepare_shiny.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
