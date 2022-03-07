#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G
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

MAINDIR="/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"

## Delete the logs and re-submit the build SPE script
cd ${CODEDIR}/04_build_spe
rm logs/build_basic_spe.txt
rm ${PROCESSEDIR}/04_build_spe/spe*.Rdata
rm logs/initial_exploration.txt
qsub 01_build_basic_spe.sh
qsub 02_initial_exploration.sh

## Run the spot QC code
cd ${CODEDIR}/07_spot_qc
rm logs/qc_metrics_and_segmentation.txt
rm logs/compare_pathology_number_to_percent.txt
rm ${PROCESSEDIR}/07_spot_qc/spe*.Rdata
qsub 01_qc_metrics_and_segmentation.sh
qsub 02_compare_pathology_number_to_percent.sh

## Run harmony and BayesSpace
cd ${CODEDIR}/08_harmony_BayesSpace
rm ${PROCESSEDIR}/08_harmony_BayesSpace/*/spe_harmony_*.rds
rm logs/preprocess_and_harmony*.txt
rm logs/BayesSpace_k_search_spe_harmony*.txt
rm logs/plot_SNN_k10.txt
rm logs/heatmaps_dlpfc_markers.txt
rm logs/fasthplus_optimal_k_*.txt
sh 01_preprocess_and_harmony.sh
sh 02_BayesSpace_k_search.sh
qsub 03_plot_SNN_k10.sh
qsub 04_heatmaps_dlpfc_markers.sh
sh 05_fasthplus_optimal_k.sh

## Run code related to comparing pathology vs BayesSpace
cd ${CODEDIR}/09_pathology_vs_BayesSpace
rm logs/bayesspace_k_pathology_boxplots.txt
rm logs/pathology_thresholds.txt
rm logs/bayesspace_pathology_barplots.txt
rm logs/label_pathology_spots.txt
qsub 01_bayesspace_k_pathology_boxplots.sh
qsub 02_pathology_thresholds.sh
qsub 03_bayesspace_pathology_barplots.sh
qsub 04_label_pathology_spots.sh

## Add future steps here
## TODO

## Delete SPE versions for shiny
rm ${CODEDIR}/0*_deploy_app*/spe*.Rdata

## Prepare the inputs for shiny (steps 05 and 06)
## (this will likely be the last step)
cd ${CODEDIR}/99_prepare_for_shiny
rm logs/prepare_shiny*.txt
qsub 01_prepare_shiny.sh

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
