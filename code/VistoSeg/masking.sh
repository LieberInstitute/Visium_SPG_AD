#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4 
#$ -N dotdotdot
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/masking_test.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/masking_test.txt
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
##echo "Sample id: VIFAD2_V10A27-106_B1"
##echo "Sample id: VIFAD3_V10T31-036_D1"
echo "****"

## load matlab
module load matlab/R2019a

## List current modules for reproducibility
module list

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/ 

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_A1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[4001,5000,7101,8100])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_A1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[4001,5000,7101,8100])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_A1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[8001,9000,7501,8500])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[8001,9000,7501,8500])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_C1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[15301,16300,9801,10800])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_D1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[9501,10500,4801,5800])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[5001,6000,9501,10500])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_B1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[5001,6000,9501,10500])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_C1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[7301,8300,9801,10800])"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_D1_segmentation.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), masking('$filename',[4001,5000,7101,8100])"


echo "**** Job ends ****"
date

