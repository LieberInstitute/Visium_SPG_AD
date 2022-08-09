#### load relevant packages ####
library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')
library('biomaRt')
library('spatialLIBD')

# library('RColorBrewer')
# library('ggplot2')
# library('fields')
# library('limma')
# library('jaffelab')
# library('janitor')
# library('lattice')
# library('org.Hs.eg.db')
# library('GenomicFeatures')

### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensemble_function.R'))

#### read in necessary input files ####
visad_enrichment_stats <- read.csv(here('code','05_deploy_app_wholegenome',
                                        'Visium_IF_AD_wholegenome_model_results_FDR5perc.csv'))
head(visad_enrichment_stats)

mostafavi_dir <- here('raw-data', 'GeneSets',
                      '1_Bulk_RNA-seq' ,'Mostafavi_et_al')
table_s3 <- read_xlsx(paste0(mostafavi_dir, '/Table_S3_M109_390_genes.xlsx'),
                      sheet = 1, col_names = TRUE, skip = 4)
#nrow
#13153

table_s3 <- table_s3 |> filter(`Module ID` == 'm109')
#390

# head(table_s3)
# # A tibble: 6 × 2
# `Module ID` `Gene Symbol`
# <chr>       <chr>
#     1 m109        RP11-742N3.1
# 2 m109        SLC39A11
# 3 m109        HEY2
# 4 m109        SMARCC1
# 5 m109        RBM4B
# 6 m109        SYNRG

table_s8 <-read_xlsx(paste0(mostafavi_dir, '/Table_S8_M109_112_genes.xlsx'),
                     sheet = 1, col_names = TRUE, skip =2 )


# # A tibble: 6 × 5
# `Gene symbol` `Degree in BN` iNs                Astrocytes Microglia
# <chr>                  <dbl> <chr>              <chr>      <chr>
#     1 CSRP1                     13 2.97               3.65       3.02
# 2 PLXNB1                    12 14.48              17.63      12.11
# 3 FAM63A                    12 2.4500000000000002 3.44       4.08
# 4 KIF1C                     11 2.99               1.67       0
# 5 CCDC85C                   11 0.34               0.71       1.44
# 6 HMG20B                    11 3.37               7.83       0.57999999999999996

table_s9 <-read_xlsx(paste0(mostafavi_dir, '/Table_S9_M109_21_genes.xlsx'),
                     sheet = 1, col_names = TRUE, skip = 4)

# A tibble: 6 × 5
# Gene   ...2  `Target Sequence`     `Broad Public ID` `Target Region`
# <chr>  <chr> <chr>                 <chr>             <chr>
#     1 BCL2L1 A     GCTCACTCTTCAGTCGGAAAT TRCN0000033499    CDS
# 2 NA     B     GTGGAACTCTATGGGAACAAT TRCN0000033500    CDS
# 3 NA     C     GTTTAGTGATGTGGAAGAGAA TRCN0000299586    CDS
# 4 NA     D     CGACGAGTTTGAACTGCGGTA TRCN0000033503    CDS
# 5 NA     E     AGAGCTTTGAACAGGATACTT TRCN0000033502    CDS
# 6 KIF5B  A     TTACAACTGTGGCCCTATTTA TRCN0000338580    3UTR




get_ensemble(table = table_s3, gene_col = "Gene Symbol")
get_ensemble(table = table_s8, gene_col = "Gene symbol")
get_ensemble(table = table_s9, gene_col = "Gene")

gene_set_enrichment(
    gene_list,
    fdr_cut = 0.1,
    modeling_results = fetch_data(type = "modeling_results"),
    model_type = names(modeling_results)[1],
    reverse = FALSE
)
#multiple ENSEMBL IDs?
