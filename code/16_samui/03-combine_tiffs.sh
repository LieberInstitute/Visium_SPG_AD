#!/bin/bash
#$ -cwd
#$ -N "combine_tiffs"
#$ -o ../../processed-data/16_samui/logs/03-combine_tiffs_$TASK_ID.log
#$ -e ../../processed-data/16_samui/logs/03-combine_tiffs_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
#$ -t 5
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load loopy/1.0.0-next.24
python 03-combine_tiffs.py

echo "**** Job ends ****"
date
