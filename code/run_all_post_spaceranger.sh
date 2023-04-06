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

MAINDIR="/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD"
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
rm logs/BayesSpace_k_search_*.txt
rm logs/plot_SNN_k10.txt
rm logs/heatmaps_dlpfc_markers.txt
rm logs/fasthplus_optimal_k_*.txt
rm logs/BayesSpace_enhance_a_split*.txt
rm logs/BayesSpace_enhance_b_compute_*.txt
rm logs/BayesSpace_enhance_c_merge*.txt
sh 01_preprocess_and_harmony.sh
sh 02_BayesSpace_k_search.sh
qsub 03_plot_SNN_k10.sh
qsub 04_heatmaps_dlpfc_markers.sh
sh 05_fasthplus_optimal_k.sh
qsub 06_BayesSpace_enhance_a_split.sh
sh 06_BayesSpace_enhance_b_compute.sh
qsub 06_BayesSpace_enhance_c_merge.sh

## Run code related to comparing pathology vs BayesSpace
cd ${CODEDIR}/09_pathology_vs_BayesSpace
rm logs/bayesspace_k_pathology_boxplots.txt
rm logs/pathology_thresholds.txt
rm logs/bayesspace_pathology_barplots.txt
rm logs/label_pathology_spots.txt
rm logs/BayesSpace_pathology_barplots_v2.txt
qsub 01_bayesspace_k_pathology_boxplots.sh
qsub 02_pathology_thresholds.sh
qsub 03_bayesspace_pathology_barplots.sh
qsub 04_label_pathology_spots.sh
qsub 05_BayesSpace_pathology_barplots_v2.sh

## Perform spatial registration
cd ${CODEDIR}/10_spatial_registration
rm logs/spatial_registration*.txt
sh 01_spatial_registration.sh

## Grey matter only analysis
cd ${CODEDIR}/11_grey_matter_only
rm logs/create_pseudobulk_data*.txt
rm logs/explore_expr_variability*.txt
rm logs/model_pathology*.txt
sh 01_create_pseudobulk_data.sh
sh 02_explore_expr_variability.sh
sh 03_model_pathology.sh

##MAGMA Analysis
cd ${CODEDIR}/12_magma/01_Jansen_2019
rm logs/pvalue_based_gene_sets.txt
rm ad_gwas_50.log
rm ad_gwas_100.log
rm ad_gwas_200.log
sh 03_pvalue_based_gene_sets.sh
sh 02_magma_step_2_AD.sh

cd ${CODEDIR}/12_magma/02_Lancet_2014
rm logs/01_magma_setup.FTD.txt
rm ftd_gwas_50.log
rm ftd_gwas_100.log
rm ftd_gwas_200.log
rm snp-wise.log
rm snp-wise.log.suppl
sh 01_magma_setup_FTD.sh
sh 02_magma_all_steps_FTD.sh

cd ${CODEDIR}/12_magma/03_Nalls_2019
rm logs/01_magma_setup.FTD.txt
rm pd_gwas_50.log
rm pd_gwas_100.log
rm pd_gwas_200.log
sh 02_magma_step_3_PD.sh

cd ${CODEDIR}/12_magma
rm logs/magma_heatmap.txt
sh magma_heatmap.sh


##TGE panel exploration
cd ${CODEDIR}/13_tge_panel


##External Gene Sets
cd ${CODEDIR}/14_external_gene_sets

## Grey matter only analysis with the Abeta micro environment
cd ${CODEDIR}/17_grey_matter_only_Abeta_microenv
rm logs/create_pseudobulk_data*.txt
rm logs/explore_expr_variability*.txt
rm logs/model_pathology*.txt
sh 01_create_pseudobulk_data.sh
sh 02_explore_expr_variability.sh
sh 03_model_pathology.sh


## Delete SPE versions for shiny
rm ${CODEDIR}/0*_deploy_app*/spe*.Rdata

## prepare to share the data through spatialLIBD
cd ${CODEDIR}/98_prepare_to_share
rm logs/prepare_to_share*.txt
qsub 01_prepare_to_share.sh

## Prepare the inputs for shiny (steps 05 and 06)
## (this will likely be the last step)
cd ${CODEDIR}/99_prepare_for_shiny
rm logs/prepare_shiny*.txt
qsub 01_prepare_shiny.sh

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
