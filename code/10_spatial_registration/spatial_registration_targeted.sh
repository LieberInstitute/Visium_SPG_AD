#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N spatial_registration_targeted
#$ -o logs/spatial_registration_targeted.$TASK_ID.txt
#$ -e logs/spatial_registration_targeted.$TASK_ID.txt
#$ -m e
#$ -t 2-28
#$ -tc 20
#$ -hold_jid BayesSpace_k_search_targeted


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: sparthib"
echo "Job id: 3115891"
echo "Job name: ${SHORT}"
echo "Hostname: compute-095.cm.cluster"
echo "Task id: 26"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_spatial_registration.R "targeted"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

