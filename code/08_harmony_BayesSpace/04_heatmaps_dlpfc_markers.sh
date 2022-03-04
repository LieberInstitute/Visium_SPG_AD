#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=6G,h_vmem=6G,h_fsize=100G
#$ -N heatmaps_dlpfc_markers
#$ -o logs/heatmaps_dlpfc_markers.txt
#$ -e logs/heatmaps_dlpfc_markers.txt
#$ -m e
#$ -hold_jid preprocess_and_harmony_wholegenome,preprocess_and_harmony_targeted
#
# Should really be BayesSpace_k_search_wholegenome,BayesSpace_k_search_targeted

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
Rscript 04_heatmaps_dlpfc_markers.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
