#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/countSpots.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/countSpots.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 2-10
#$ -tc 5

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/listOFfiles.txt | awk "NR==${SGE_TASK_ID}")"
echo "****"

## List current modules for reproducibility
module list

mask=$(awk 'BEGIN {FS="\t"} {print $1}' /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/listOFfiles.txt | awk "NR==${SGE_TASK_ID}")
jsonname=$(awk 'BEGIN {FS="\t"} {print $2}' /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/listOFfiles.txt | awk "NR==${SGE_TASK_ID}")
posname=$(awk 'BEGIN {FS="\t"} {print $3}' /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/listOFfiles.txt | awk "NR==${SGE_TASK_ID}")

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')" 

echo "**** Job ends ****"
date

