#!/bin/bash
#$ -cwd
#$ -l mem_free=60G,h_vmem=60G,h_fsize=100G
#$ -N gene_pairs_path
#$ -o logs/01_gene_pairs_path.$TASK_ID.txt
#$ -e logs/01_gene_pairs_path.$TASK_ID.txt
#$ -m e
#$ -t 1-7
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_gene_pairs_path.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
