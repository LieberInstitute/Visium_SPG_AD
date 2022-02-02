#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -pe local 4
#$ -N preprocess_and_harmony_spe_targeted_postqc
#$ -o logs/preprocess_and_harmony_spe_targeted_postqc.txt
#$ -e logs/preprocess_and_harmony_spe_targeted_postqc.txt
#$ -m e
#$ -hold_jid qc_metrics_and_segmentation

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
Rscript 01_preprocess_and_harmony.R -s spe_targeted_postqc

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


