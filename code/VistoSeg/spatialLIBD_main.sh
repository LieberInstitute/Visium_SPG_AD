#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/spatialLIBD.txt
#$ -e /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/logs/spatialLIBD.txt
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
## List current modules for reproducibility
module list

matlab -nodesktop -nosplash -nojvm -r "spatialLIBD_main"

echo "**** Job ends ****"
date
