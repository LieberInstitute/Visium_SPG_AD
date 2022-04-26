#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -N explore_expr_variability_wholegenome
#$ -o logs/explore_expr_variability_wholegenome.txt
#$ -e logs/explore_expr_variability_wholegenome.txt
#$ -m e
#$ -hold_jid create_pseudobulk_data_wholegenome

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
Rscript 02_explore_expr_variability.R -s wholegenome

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


