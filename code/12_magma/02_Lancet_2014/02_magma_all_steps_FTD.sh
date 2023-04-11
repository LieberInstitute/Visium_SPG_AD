#!/bin/bash
#$ -cwd
#$ -N magma_all_steps_FTD
#$ -o ./logs/magma_all_steps_FTD.o
#$ -e ./logs/magma_all_steps_FTD.e
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
SNPLOC=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/magma/02_Lancet_2014/FTD_Lancet2014.snploc
STEP1OUT=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/magma/02_Lancet_2014/FTD_Lancet2014_ITC

BFILE=/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
ANNOT=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/magma/02_Lancet_2014/FTD_Lancet2014_ITC.genes.annot

setcol=1
genecol=2

here=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/12_magma
gs_200=$here/pvalues_top_200.txt
gs_50=$here/pvalues_top_50.txt
gs_100=$here/pvalues_top_100.txt
fdr_set=$here/fdr_gene_set.txt
SUMMSTATS=$here/02_Lancet_2014/FTD-IFGC-and-rsID-ADDED.tab
STEP3OUT_200=$here/02_Lancet_2014/step3_results/ftd_gwas_200
STEP3OUT_50=$here/02_Lancet_2014/step3_results/ftd_gwas_50
STEP3OUT_100=$here/02_Lancet_2014/step3_results/ftd_gwas_100
STEP3OUT_FDR=$here/02_Lancet_2014/step3_results/ftd_gwas_fdr

## Step 1 - Annotation (SNP : gene mapping)
magma --annotate window=35,10 --snp-loc $SNPLOC --gene-loc $ANNO --out $STEP1OUT

## Step 2 - Gene analysis (from SNP-level summary stats)
magma --bfile $BFILE --gene-annot $ANNOT --pval $SUMMSTATS use=rsID,pValue ncol=N_effective --gene-model ${model} --out /dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/code/magma/02_Lancet_2014/step2_results/${model}


## Step 3 - Gene set analyses (using gene-level output) - Self-contained
#for top 200
magma --gene-results $here/02_Lancet_2014/snp-wise.genes.raw --set-annot $gs_200 gene-col=${genecol} set-col=${setcol} --out $STEP3OUT_200
#for top 100
magma --gene-results $here/02_Lancet_2014/snp-wise.genes.raw --set-annot $gs_100 gene-col=${genecol} set-col=${setcol} --out $STEP3OUT_100
#for top 50
magma --gene-results $here/02_Lancet_2014/snp-wise.genes.raw --set-annot $gs_50 gene-col=${genecol} set-col=${setcol} --out $STEP3OUT_50

#fdr geneset
magma --gene-results $here/02_Lancet_2014/snp-wise.genes.raw --set-annot $fdr_set gene-col=${genecol} set-col=${setcol} --out $STEP3OUT_FDR


echo "**** Job ends ****"
date
