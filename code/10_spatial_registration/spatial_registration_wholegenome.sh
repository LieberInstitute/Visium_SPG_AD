#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N spatial_registration_wholegenome
#$ -o logs/spatial_registration_wholegenome.$TASK_ID.txt
#$ -e logs/spatial_registration_wholegenome.$TASK_ID.txt
#$ -m e
#$ -t 2-28
#$ -tc 20


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: lcollado"
echo "Job id: "
echo "Job name: "
echo "Hostname: compute-122.cm.cluster"
echo "Task id: "

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_spatial_registration.R -s wholegenome

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

