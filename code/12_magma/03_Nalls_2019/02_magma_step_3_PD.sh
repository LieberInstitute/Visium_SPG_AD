#!/bin/bash
#$ -cwd
#$ -N magma_step3_PD
#$ -o ./logs/magma_step3_PD.o
#$ -e ./logs/magma_step3_PD.e
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
# ANNO=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/annotation/GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc
# BFILE=/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
#
setcol=1
genecol=2

here=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/12_magma
gs_200=$here/pvalues_top_200.txt
gs_50=$here/pvalues_top_50.txt
gs_100=$here/pvalues_top_100.txt
fdr_set=$here/fdr_gene_set.txt
#SUMMSTATS=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/GWAS_Results/nallsEtAl2019_excluding23andMe_allVariants_Neff-and-rsID-ADDED-MNT.tab


##PD Step 3

## Step 3 - Gene set analyses (using gene-level output)

#top 200
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_200 gene-col=${genecol} set-col=${setcol} --out $here/03_Nalls_2019/pd_gwas_200

#top 100
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_100 gene-col=${genecol} set-col=${setcol} --out $here/03_Nalls_2019/pd_gwas_100

#top 50
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_50 gene-col=${genecol} set-col=${setcol} --out $here/03_Nalls_2019/pd_gwas_50

#fdr geneset
magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $fdr_set gene-col=${genecol} set-col=${setcol} --out $here/03_Nalls_2019/pd_gwas_fdr



