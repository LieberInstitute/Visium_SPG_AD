#!/bin/bash
#$ -cwd
#$ -N "SPG_images"
#$ -o ../../processed-data/16_samui/logs/01-SPG_images_$TASK_ID.log
#$ -e ../../processed-data/16_samui/logs/01-SPG_images_$TASK_ID.log
#$ -l mf=10G,h_vmem=10G
#$ -t 1
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load loopy/1.0.0-next.15
python 01-SPG_images.py

echo "**** Job ends ****"
date
