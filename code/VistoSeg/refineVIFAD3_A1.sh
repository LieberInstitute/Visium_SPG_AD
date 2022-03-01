#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/refineVIFAD3_A1.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/refineVIFAD3_A1.txt
#$ -m e
#$ -M madhavitippani28@gmail.com

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1.mat"
echo "****"

## List current modules for reproducibility
module list
scp /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1.mat /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/VIFAD3_A1/
fname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/VIFAD3_A1/VIFAD3_V10T31-036_A1.mat'
toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/'
matlab -nodesktop -nosplash -r "addpath(genpath('$toolbox')), refineIF('$fname',[13501,14500,8701,9700],[],[],0.15,0.1)"

echo "**** Job ends ****"
date

