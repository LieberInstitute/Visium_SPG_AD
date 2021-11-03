#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=200G
#$ -N update_backup_Visium_IF_AD
#$ -o logs/update_backup.txt
#$ -e logs/update_backup.txt
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

## Set restrictive permissions
umask 077

## Update the backup location
rsync -avh /dcs04/lieber/lcolladotor/with10x_LIBD001 /dcl02/lieber/lcolladotor/backup_with10x_LIBD001

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/