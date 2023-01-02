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

gs_200=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/12_magma/pvalues_top_200.txt
gs_50=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/12_magma/pvalues_top_50.txt
gs_100=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/12_magma/pvalues_top_100.txt
fdr_set=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/12_magma/fdr_gene_set.txt
SUMMSTATS=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/12_magma/GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt


##AD
## Step 1 - Annotation (SNP : gene mapping)
# magma --annotate window=35,10 --snp-loc ./GWAS_Results/Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc --gene-loc $ANNO --out SNP_Data/AD_Jansen2019_LC

## Step 2 - Gene analysis (from SNP-level summary stats)
# magma --bfile $BFILE --gene-annot SNP_Data/AD_Jansen2019_LC.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/AD_Jansen2019_LC_${model}


## Step 3 - Gene set analyses (using gene-level output)

#top 200
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.raw --set-annot $gs_200 gene-col=${genecol} set-col=${setcol} --out code/12_magma/01_Jansen_2019/ad_gwas_200

#top 100
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.raw --set-annot $gs_100 gene-col=${genecol} set-col=${setcol} --out code/12_magma/01_Jansen_2019/ad_gwas_100

#top 50
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.raw --set-annot $gs_50 gene-col=${genecol} set-col=${setcol} --out code/12_magma/01_Jansen_2019/ad_gwas_50

#fdr geneset
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.raw --set-annot $fdr_set gene-col=${genecol} set-col=${setcol} --out code/12_magma/01_Jansen_2019/ad_gwas_fdr

echo "**** Job ends ****"
date
