#!/bin/bash
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/InFormStitch.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/logs/InFormStitch.txt
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

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/VistoSeg/'

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD1_V10A27-004/20210204_Visium-IF_Scan1_*_component_data.tif'
fname='VIFAD1_V10A27-004'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Abeta'; O{3} = 'pTau'; O{4} = 'GFAP'; O{5} = 'MAP2'; O{6} = 'Lipofuscin'; InFormStitch('$filename',O,4,'$fname')"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD2_V10A27-106/20210331_Visium-IF2nd_Scan1_*_component_data.tif'
fname='VIFAD2_V10A27-106'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Abeta'; O{3} = 'pTau'; O{4} = 'GFAP'; O{5} = 'MAP2'; O{6} = 'Lipofuscin'; InFormStitch('$filename',O,4,'$fname')"

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/InForm/VIFAD3_V10T31-036/20210405_VIFAD_3rd_Scan2_*_component_data.tif'
fname='VIFAD3_V10T31-036'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Abeta'; O{3} = 'pTau'; O{4} = 'GFAP'; O{5} = 'MAP2'; O{6} = 'Lipofuscin'; InFormStitch('$filename',O,5,'$fname')"

echo "**** Job ends ****"
date

