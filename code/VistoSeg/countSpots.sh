#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/countSpots_8.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/countSpots_8.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 12
#$ -tc 3

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/'

mask='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_D1_segmentation.mat'
jsonname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/spaceranger/V10A27106_D1_Br3880/outs/spatial/scalefactors_json.json'
posname='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/spaceranger/V10A27106_D1_Br3880/outs/spatial/tissue_positions_list.csv'

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')" 

echo "**** Job ends ****"
date

