#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4 
#$ -N dotdotdot
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/dotdotdot_test.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/dotdotdot_test.txt
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

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/'

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_A1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[8301,9300,9501,10500],0.3)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_D1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[4001,5000,7101,8100],0.4)"
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[4001,5000,11501,12500],0.4)"
	
##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_A1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[12501,13500,5801,6800],0.2,0)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[8701,9700,6801,7800],0.3,0)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_C1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[15301,16300,9801,10800],0.4)"
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[15001,16000,7601,8600],0.4)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_D1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[4801,5800,8001,9000],0.3,0.09)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[5001,6000,9501,10500],0.7)"
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[13501,14500,8701,9700],0.7)"

##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_B1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[5001,6000,9501,10500],0.7)"
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[7001,8000,9501,10500],0.7)"


##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_C1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[7301,8300,9801,10800],0.7)"
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[13001,14000,10201,11200],0.7)"


##filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_D1.mat'
##matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), rnascope_human('$filename',[4001,5000,7101,8100],0.7)"
##[4251,5250,7501,8500]


echo "**** Job ends ****"
date

