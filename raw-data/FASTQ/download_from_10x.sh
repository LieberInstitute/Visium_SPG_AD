#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=200G
#$ -N download_from_10x
#$ -o logs/download_from_10x.txt
#$ -e logs/download_from_10x.txt
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
module load conda_R

## List current modules for reproducibility
module list

## Download data and uncompress
wget https://ty.10xgenomics.com/6ef33f518641b336/Leiber_Fastqs.tar
tar -xvf Leiber_Fastqs.tar
mv Leiber_Fastqs FASTQ_10x_2021-05-03

## Hm... uncompress the compressed files
tar -xvf FASTQ_10x_2021-05-03/Leiber/fastqs.tar.gz

## Re-organize and delete files we don't need
mv FASTQ_10x_2021-05-03/Leiber/mnt/analysis/marsoc/pipestances/H5GKNDSX2/BCL_PROCESSOR_PD/H5GKNDSX2/2020.0623.2-0/outs/fastq_path/* FASTQ_10x_2021-05-03/
rm -fr FASTQ_10x_2021-05-03/Leiber
rm Leiber_Fastqs.tar

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
