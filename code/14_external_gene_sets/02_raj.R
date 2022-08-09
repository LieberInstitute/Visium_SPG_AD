#### load relevant packages ####
library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')



# Table S2: split by "Trait". There's 3 of them.
# Direction available: Beta, but we won't use it since there can be an up and a down intron for the same gene.
# Statistics available: FDR
# Table S3:
#     Direction available: None
# Statistics available: Both FDR (P-value_Benjamini-Hochberg) and Bonferroni (P-value_BonferroniAdjusted) adjusted.
# Could subset each to < 0.05 (so 2 sets total).


### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensemble_function.R'))


#### read in necessary input files ####
table_s2 <- read_excel("raw-data/GeneSets/1_Bulk_RNA-seq/Raj et al/Table S2.xlsx")
head(table_s2)
# intronic_cluster_id            cluster     chr    start      end gene_id    Beta     SE `Z-score`   `P-value`    FDR Trait
# <chr>                          <chr>     <dbl>    <dbl>    <dbl> <chr>     <dbl>  <dbl>     <dbl>       <dbl>  <dbl> <chr>
#     1 10_3146136_3147307_clu_8247    clu_8247     10  3146136  3147307 PFKP    -0.0793 0.0160     -4.97 0.000000967 0.0158 NEURITIC PLAQUES
# 2 10_3146980_3147307_clu_8247    clu_8247     10  3146980  3147307 PFKP     0.119  0.0231      5.14 0.000000416 0.0134 NEURITIC PLAQUES
# 3 10_3147351_3147585_clu_8247    clu_8247     10  3147351  3147585 PFKP     0.0898 0.0175      5.13 0.000000430 0.0134 NEURITIC PLAQUES


table_s3 <- read_excel("raw-data/GeneSets/1_Bulk_RNA-seq/Raj et al/Table S3.xlsx")
# intronic_cluster gene_id       `P-value` `P-value_Benjamini-Hochberg` `P-value_BonferroniAdjusted`
# <chr>            <chr>             <dbl>                        <dbl>                        <dbl>
#     1 chr10:clu_8247   PFKP           1.79e-28                     4.95e-24                     4.97e-24
# 2 chr14:clu_18324  NDRG2          2.03e-23                     2.81e-19                     5.62e-19
# 3 chr19:clu_21882  CTD-2527I21.4  2.74e-20                     2.53e-16                     7.59e-16
# 4 chr7:clu_28844   BCL7B          6.95e-19                     4.81e-15                     1.92e-14


