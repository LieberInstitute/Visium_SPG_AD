#!/bin/bash
#$ -cwd
#$ -N "samui_link"
#$ -o ../../processed-data/16_samui/logs/04-samui_link.log
#$ -e ../../processed-data/16_samui/logs/04-samui_link.log
#$ -l mf=3G,h_vmem=3G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 04-samui_link.R

echo "**** Job ends ****"
date
