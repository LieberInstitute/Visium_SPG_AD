

echo "**** Job starts ****"
date

## List current modules for reproducibility
module list

## Load MAGMA
module load magma/1.10

#magma --gene-results /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_ad gene-col=${genecol} set-col=${setcol} --out code/magma/03_Nalls_2019/pd_gwas_200
