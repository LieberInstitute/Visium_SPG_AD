#!/bin/bash
#$ -cwd
#$ -N magma_all_steps_AD
#$ -o ./logs/magma_all_steps_AD.o
#$ -e ./logs/magma_all_steps_AD.e
#$ -l bluejay,mem_free=16G,h_vmem=20G

#borrowed from Matt's script here https://github.com/lmweber/locus-c/blob/main/
#code/magma/magma-gsa_allSteps_AD_LCclusterMarkers.sh
echo "**** Job starts ****"
date

## List current modules for reproducibility
module list

## Load MAGMA
module load magma/1.10

## Set some variables/paths
model="snp-wise"
ANNO=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/annotation/GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc
BFILE=/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur

setcol=1
genecol=2

gs_ad=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/magma/pvalues_top_200.txt
SUMMSTATS=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt
##

## Step 1 - Annotation (SNP : gene mapping)
# magma --annotate window=35,10 --snp-loc ./GWAS_Results/Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc --gene-loc $ANNO --out SNP_Data/AD_Jansen2019_LC

## Step 2 - Gene analysis (from SNP-level summary stats)
# magma --bfile $BFILE --gene-annot SNP_Data/AD_Jansen2019_LC.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/AD_Jansen2019_LC_${model}


## Step 3 - Gene set analyses (using gene-level output)
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.raw --set-annot $gs_ad gene-col=${genecol} set-col=${setcol} --out code/magma/01_Jansen_2019/ad_gwas_200
#magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_ad gene-col=${genecol} set-col=${setcol} --out code/magma/03_Nalls_2019/pd_gwas_200


echo "**** Job ends ****"
date
