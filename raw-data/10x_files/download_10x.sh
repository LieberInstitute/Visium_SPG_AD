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
## for less verbose output: https://www.davekb.com/browse_computer_tips:wget_disable_progress_bar:txt
wget -nv https://ty.10xgenomics.com/91c621c7b18d8b86/Lieber_Transfer.tar
wget -nv https://ty.10xgenomics.com/42e0d8f4a79cea20/Lieber_Transfer_10x_Alignments.tar

## Uncompress
tar -xvf Lieber_Transfer.tar
tar -xvf Lieber_Transfer_10x_Alignments.tar

## Clean up
rm Lieber_Transfer.tar
rm Lieber_Transfer_10x_Alignments.tar

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
