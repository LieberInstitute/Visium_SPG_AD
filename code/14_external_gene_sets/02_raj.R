#### load relevant packages ####
# library("sgejobs")

# sgejobs::job_single(
#     "raj",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 02_raj.R"
# )

library("readxl")
library("spatialLIBD")
library("dplyr")
library("sessioninfo")
library("here")
library("scran")
library("purrr")



# Table S2: split by "Trait". There's 3 of them.
# Statistics available: FDR
# Table S3:
#     Direction available: None
# Statistics available: Both FDR (P-value_Benjamini-Hochberg) and Bonferroni (P-value_BonferroniAdjusted) adjusted.
# Could subset each to < 0.05 (so 2 sets total).


### load get_ensemble function
source(here("code/14_external_gene_sets/get_ensembl_function.R"))

#### load modeling results ####
load(here(
    "processed-data", "11_grey_matter_only", "wholegenome",
    "Visium_SPG_AD_modeling_results.Rdata"
))

#### read in external gene sets  ####
table_s2 <- read_excel(here("raw-data", "GeneSets", "1_Bulk_RNA-seq", "Raj et al", "Table S2.xlsx"))

head(table_s2)
# intronic_cluster_id            cluster     chr    start      end gene_id    Beta     SE `Z-score`   `P-value`    FDR Trait
# <chr>                          <chr>     <dbl>    <dbl>    <dbl> <chr>     <dbl>  <dbl>     <dbl>       <dbl>  <dbl> <chr>
#     1 10_3146136_3147307_clu_8247    clu_8247     10  3146136  3147307 PFKP    -0.0793 0.0160     -4.97 0.000000967 0.0158 NEURITIC PLAQUES
# 2 10_3146980_3147307_clu_8247    clu_8247     10  3146980  3147307 PFKP     0.119  0.0231      5.14 0.000000416 0.0134 NEURITIC PLAQUES
# 3 10_3147351_3147585_clu_8247    clu_8247     10  3147351  3147585 PFKP     0.0898 0.0175      5.13 0.000000430 0.0134 NEURITIC PLAQUES

unique(table_s2$Trait)
# "NEURITIC PLAQUES",  "AMYLOID", "Tangles"


table_s2 <- table_s2 |> dplyr::filter(FDR < 0.1)
unique(table_s2 |> dplyr::filter(!is.na(gene_ensembl_id)) |> dplyr::select(gene_id))
table_s2_np <- table_s2 |> dplyr::filter(Trait == "NEURITIC PLAQUES")
table_s2_am <- table_s2 |> dplyr::filter(Trait == "AMYLOID")
genes_s2_am <- get_ensembl(table_s2_am, gene_id, "gene_id")
table_s2_ta <- table_s2 |> dplyr::filter(Trait == "Tangles")
genes_s2_ta <- get_ensembl(table_s2_ta, gene_id, "gene_id")


# gene_id
# 1     AC015936
# 2     AC018730
# 5     AC068057
# 8     AC174470
# 9        ATP5J
# 10    C19orf26
# 11      FAM73B
# 12    FLJ27365
# 15   KIAA1211L
# 16       MRP63
# 17        PIDD
# 18  RP11-123K3
# 19 RP11-463D19
# 20 RP11-637O19
# 21  RP11-85M11
# 23   RP6-109B7
# 24        SARS
# 25      SUZ12P

nrow(table_s2)
# [1] 234

# > nrow(table_s2_np)
# [1] 2
# > nrow(table_s2_am)
# [1] 65
# > nrow(table_s2_ta)
# [1] 167

# load table 3
table_s3 <- read_excel(here("raw-data", "GeneSets", "1_Bulk_RNA-seq", "Raj et al", "Table S3.xlsx"))
# intronic_cluster      gene_id       `P-value` `P-value_Benjamini-Hochberg` `P-value_BonferroniAdjusted`
# <chr>            <chr>             <dbl>                        <dbl>                        <dbl>
#     1 chr10:clu_8247   PFKP           1.79e-28                     4.95e-24                     4.97e-24
# 2 chr14:clu_18324  NDRG2          2.03e-23                     2.81e-19                     5.62e-19
# 3 chr19:clu_21882  CTD-2527I21.4  2.74e-20                     2.53e-16                     7.59e-16
# 4 chr7:clu_28844   BCL7B          6.95e-19                     4.81e-15                     1.92e-14

table_s3 <- table_s3 |> dplyr::filter(`P-value_BonferroniAdjusted` < 0.1)
# > nrow(table_s3)
# [1] 99

table_s3 <- get_ensembl(table_s3, gene_id, "gene_id")

# > nrow(table_s3 |> dplyr::filter(is.na(gene_ensembl_id)))
# [1] 4


raj_geneList <- list(
    raj_table_2_am = genes_s2_am,
    raj_table_2_ta = genes_s2_ta,
    raj_table_3 = table_s3
)






#### calculate enrichment #####
raj_enrichment <- gene_set_enrichment(
    raj_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

raj_depleted <- gene_set_enrichment(
    raj_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment", reverse = TRUE
)


##### Enrichment plotting #####
# dir.create(here("plots", "14_external_gene_sets"))
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/02_raj_enriched.pdf"), width = 11)

gene_set_enrichment_plot(
    raj_enrichment,
    xlabs = unique(raj_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(raj_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

pdf(paste0(output_dir, "/02_raj_depleted.pdf"), width = 11)

gene_set_enrichment_plot(
    raj_depleted,
    xlabs = unique(raj_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(raj_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
