#!/bin/bash
#$ -cwd
#$ -N "spe_to_anndata"
#$ -o ../../processed-data/16_samui/logs/02-spe_to_anndata.log
#$ -e ../../processed-data/16_samui/logs/02-spe_to_anndata.log
#$ -l mf=70G,h_vmem=70G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 02-spe_to_anndata.R

echo "**** Job ends ****"
date
