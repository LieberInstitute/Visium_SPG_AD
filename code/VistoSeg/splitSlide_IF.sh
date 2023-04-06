#!/bin/bash
#$ -cwd;
#$ -pe local 10
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/splitSlide.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/splitSlide.txt
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

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/'

fname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD1_V10A27-004/VIFAD1_V10A27-004.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fnamefunction ')" 
mv /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD1_V10A27-004/VIFAD1_V10A27-004_*1* /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas

fname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD2_V10A27-106/VIFAD2_V10A27-106.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname')"
mv /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD2_V10A27-106/VIFAD2_V10A27-106_*1* /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas

fname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD3_V10T31-036/VIFAD3_V10T31-036.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname')"
mv /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD3_V10T31-036/VIFAD3_V10T31-036_*1* /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas

echo "**** Job ends ****"
date

