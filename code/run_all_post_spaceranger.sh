#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -N VisiumIFAD_run_all_post_spaceranger
#$ -o logs/run_all_post_spaceranger.txt
#$ -e logs/run_all_post_spaceranger.txt
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

CODEDIR="/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code"

## Delete the logs and re-submit the build SPE script
cd ${CODEDIR}/04_build_spe
rm logs/build_basic_spe.txt
qsub build_basic_spe.sh

## Delete SPE versions for shiny
rm ${CODEDIR}/0*_deploy_app*/spe*.Rdata

## Run the spot QC code
cd ${CODEDIR}/07_spot_qc
rm logs/qc_metrics_and_segmentation.txt
rm logs/subset_data_shiny.txt
rm logs/compare_pathology_number_to_percent.txt
qsub qc_metrics_and_segmentation.sh
qsub subset_data_shiny.sh
qsub compare_pathology_number_to_percent.sh

## Run harmony and BayesSpace
# cd ${CODEDIR}/08_harmony_BayesSpace
# rm logs/preprocess_and_harmony*.txt
# rm logs/BayesSpace_k_search_spe_harmony*.txt
# rm logs/plot_SNN_k10.txt
# sh preprocess_and_harmony.sh
# sh BayesSpace_k_search.sh
# qsub plot_SNN_k10.sh

## Add future steps here


echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
