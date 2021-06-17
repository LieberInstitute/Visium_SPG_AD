#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=200G
#$ -N download_10x
#$ -o logs/download_10x.txt
#$ -e logs/download_10x.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## Download the data
wget https://ty.10xgenomics.com/083ad7b0088d8912/Lieber_Transfer.tar
wget https://ty.10xgenomics.com/faf16fa65919ce60/Lieber_Transfer_10x_Alignments.tar

## Uncompress
tar -xvf Lieber_Transfer.tar
tar -xvf Lieber_Transfer_10x_Alignments.tar

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
