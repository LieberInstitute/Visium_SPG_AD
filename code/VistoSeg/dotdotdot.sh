#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -pe local 4 
#$ -N dotdotdot
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/dotdotdot_VIFAD2_B1.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/dotdotdot_VIFAD2_B1.txt
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
echo "Sample id: VIFAD2_V10A27-106_B1"
echo "****"

## load matlab
module load matlab/R2019a

## List current modules for reproducibility
module list

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/'

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[8001,9000,7501,8000])"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename')"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_C1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename')"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_D1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename')"


echo "**** Job ends ****"
date

