library("sgejobs")
# sgejobs::job_single(
#     "mostafavi",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 01_mostafavi.R",
#     create_logdir = TRUE
# )


#### load relevant packages ####

library("readxl")
library("sessioninfo")
library("here")
library("spatialLIBD")
library("dplyr")
library("scran")
library("purrr")


### load get_ensemble function
source(here("code/14_external_gene_sets/get_ensembl_function.R"))

#### read in necessary input files ####
load(here(
    "processed-data", "11_grey_matter_only", "wholegenome",
    "Visium_IF_AD_modeling_results.Rdata"
))

mostafavi_dir <- here(
    "raw-data", "GeneSets",
    "1_Bulk_RNA-seq", "Mostafavi_et_al"
)


table_s3 <- read_xlsx(paste0(mostafavi_dir, "/Table_S3_M109_390_genes.xlsx"),
    sheet = 1, col_names = TRUE, skip = 4
)
# nrow
# 13153

table_s3 <- table_s3 |> dplyr::filter(`Module ID` == "m109")
# 390

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

table_s8 <- read_xlsx(paste0(mostafavi_dir, "/Table_S8_M109_112_genes.xlsx"),
    sheet = 1, col_names = TRUE, skip = 2
)
# A tibble: 112 × 5
# `Gene symbol` `Degree in BN` iNs                Astrocytes Microglia
# <chr>                  <dbl> <chr>              <chr>      <chr>
#     1 CSRP1                     13 2.97               3.65       3.02
# 2 PLXNB1                    12 14.48              17.63      12.11
# 3 FAM63A                    12 2.4500000000000002 3.44       4.08
# 4 KIF1C                     11 2.99               1.67       0
# 5 CCDC85C                   11 0.34               0.71       1.44
# 6 HMG20B                    11 3.37               7.83       0.57999999999999996

table_s9 <- read_xlsx(paste0(mostafavi_dir, "/Table_S9_M109_21_genes.xlsx"),
    sheet = 1, col_names = TRUE, skip = 4
)

table_s9 <- table_s9 |> na.omit(Gene)
nrow(table_s9)

# A tibble: 21 × 5
# Gene   ...2  `Target Sequen…` `Broad Public …` `Target Region`
# <chr>  <chr> <chr>            <chr>            <chr>
#     1 BCL2L1 A     GCTCACTCTTCAGTC… TRCN0000033499   CDS
# 2 KIF5B  A     TTACAACTGTGGCCC… TRCN0000338580   3UTR
# 3 ACIN1  A     AGCAAGATGAGCTGG… TRCN0000236236   CDS
# 4 ITPK1  A     CGAGATGGCTATCGT… TRCN0000037695   CDS
# 5 KANSL1 A     GTTCATCCTGTTCTA… TRCN0000369821   CDS
# 6 INPPL1 A     CCACCCAAGAACAGC… TRCN0000311650   CDS


#### create gene lists # ####
table_s3_genes <- get_ensembl(table = table_s3, `Gene Symbol`, "Gene Symbol")
table_s8_genes <- get_ensembl(table = table_s8, `Gene symbol`, "Gene symbol")
table_s9_genes <- get_ensembl(table = table_s9, Gene, "Gene")




mostafavi_geneList <- list(
    mostafavi_table_3 = table_s3_genes,
    mostafavi_table_8 = table_s8_genes,
    mostafavi_table_9 = table_s9_genes
)

#### calculate enrichment #####
mostafavi_enrichment <- gene_set_enrichment(
    mostafavi_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment", reverse = FALSE
)


mostafavi_depleted <- gene_set_enrichment(
    mostafavi_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)

#> mostafavi_enrichment
# OR       Pval      test                ID model_type fdr_cut
# 1   0.000000 1.00000000      none mostafavi_table_3 enrichment     0.1
# 2   0.000000 1.00000000      none mostafavi_table_8 enrichment     0.1
# 3   0.000000 1.00000000      none mostafavi_table_9 enrichment     0.1
# 4   5.278728 0.19062682       Ab+ mostafavi_table_3 enrichment     0.1
# 5  18.018732 0.06119092       Ab+ mostafavi_table_8 enrichment     0.1
# 6   0.000000 1.00000000       Ab+ mostafavi_table_9 enrichment     0.1
# 7   0.000000 0.63325367  next_Ab+ mostafavi_table_3 enrichment     0.1
# 8   0.000000 1.00000000  next_Ab+ mostafavi_table_8 enrichment     0.1
# 9   0.000000 1.00000000  next_Ab+ mostafavi_table_9 enrichment     0.1
# 10  0.000000 1.00000000       pT+ mostafavi_table_3 enrichment     0.1
# 11  0.000000 1.00000000       pT+ mostafavi_table_8 enrichment     0.1
# 12  0.000000 1.00000000       pT+ mostafavi_table_9 enrichment     0.1
# 13  0.000000 1.00000000  next_pT+ mostafavi_table_3 enrichment     0.1
# 14  0.000000 1.00000000  next_pT+ mostafavi_table_8 enrichment     0.1
# 15  0.000000 1.00000000  next_pT+ mostafavi_table_9 enrichment     0.1
# 16  0.000000 1.00000000      both mostafavi_table_3 enrichment     0.1
# 17  0.000000 1.00000000      both mostafavi_table_8 enrichment     0.1
# 18  0.000000 1.00000000      both mostafavi_table_9 enrichment     0.1
# 19  0.000000 1.00000000 next_both mostafavi_table_3 enrichment     0.1
# 20  0.000000 1.00000000 next_both mostafavi_table_8 enrichment     0.1
# 21  0.000000 1.00000000 next_both mostafavi_table_9 enrichment     0.1

# > mostafavi_depleted
# OR        Pval      test                ID model_type fdr_cut
# 1  0.0000000 1.000000000      none mostafavi_table_3  depletion     0.1
# 2  0.0000000 1.000000000      none mostafavi_table_8  depletion     0.1
# 3  0.0000000 1.000000000      none mostafavi_table_9  depletion     0.1
# 4  1.1454394 0.595813046       Ab+ mostafavi_table_3  depletion     0.1
# 5  1.0456800 0.812156912       Ab+ mostafavi_table_8  depletion     0.1
# 6  0.0000000 0.615843932       Ab+ mostafavi_table_9  depletion     0.1
# 7  0.6046274 0.005802193  next_Ab+ mostafavi_table_3  depletion     0.1
# 8  0.7404868 0.383057105  next_Ab+ mostafavi_table_8  depletion     0.1
# 9  0.0000000 0.045245738  next_Ab+ mostafavi_table_9  depletion     0.1
# 10 0.0000000 1.000000000       pT+ mostafavi_table_3  depletion     0.1
# 11 0.0000000 1.000000000       pT+ mostafavi_table_8  depletion     0.1
# 12 0.0000000 1.000000000       pT+ mostafavi_table_9  depletion     0.1
# 13 0.0000000 1.000000000  next_pT+ mostafavi_table_3  depletion     0.1
# 14 0.0000000 1.000000000  next_pT+ mostafavi_table_8  depletion     0.1
# 15 0.0000000 1.000000000  next_pT+ mostafavi_table_9  depletion     0.1
# 16 0.0000000 1.000000000      both mostafavi_table_3  depletion     0.1
# 17 0.0000000 1.000000000      both mostafavi_table_8  depletion     0.1
# 18 0.0000000 1.000000000      both mostafavi_table_9  depletion     0.1
# 19 0.0000000 1.000000000 next_both mostafavi_table_3  depletion     0.1
# 20 0.0000000 1.000000000 next_both mostafavi_table_8  depletion     0.1
# 21 0.0000000 1.000000000 next_both mostafavi_table_9  depletion     0.1


##### Enrichment plotting #####
# dir.create(here("plots", "14_external_gene_sets"))
output_dir <- here("plots", "14_external_gene_sets")


pdf(paste0(output_dir, "/01_mostafavi_enriched.pdf"), width = 9)

gene_set_enrichment_plot(
    mostafavi_enrichment,
    xlabs = unique(mostafavi_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(mostafavi_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

pdf(paste0(output_dir, "/01_mostafavi_depleted.pdf"), width = 9)

gene_set_enrichment_plot(
    mostafavi_depleted,
    xlabs = unique(mostafavi_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(mostafavi_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()





# > sessionInfo()
# R version 4.2.0 beta (2022-04-11 r82151)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.3
#
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#
# locale:
#     [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# attached base packages:
#     [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#     [1] spatialLIBD_1.7.15          SpatialExperiment_1.5.4     scran_1.23.1                scuttle_1.5.1               SingleCellExperiment_1.17.2
# [6] SummarizedExperiment_1.25.3 Biobase_2.55.2              GenomicRanges_1.47.6        GenomeInfoDb_1.31.7         IRanges_2.29.1
# [11] S4Vectors_0.33.17           BiocGenerics_0.41.2         MatrixGenerics_1.7.0        matrixStats_0.61.0          here_1.0.1
# [16] sessioninfo_1.2.2           dplyr_1.0.8                 readxl_1.4.0                biomaRt_2.52.0
#
# loaded via a namespace (and not attached):
#     [1] utf8_1.2.2                    R.utils_2.11.0                tidyselect_1.1.2              colorout_1.2-2                RSQLite_2.2.12
# [6] AnnotationDbi_1.57.1          htmlwidgets_1.5.4             grid_4.2.0                    BiocParallel_1.29.20          DropletUtils_1.15.2
# [11] munsell_0.5.0                 ScaledMatrix_1.3.0            codetools_0.2-18              statmod_1.4.36                DT_0.22
# [16] withr_2.5.0                   colorspace_2.0-3              filelock_1.0.2                config_0.3.1                  knitr_1.38
# [21] rstudioapi_0.13               shinyWidgets_0.6.4            GenomeInfoDbData_1.2.7        bit64_4.0.5                   rhdf5_2.39.6
# [26] rprojroot_2.0.3               vctrs_0.4.1                   generics_0.1.2                xfun_0.30                     BiocFileCache_2.3.4
# [31] R6_2.5.1                      doParallel_1.0.17             ggbeeswarm_0.6.0              rsvd_1.0.5                    locfit_1.5-9.5
# [36] fields_13.3                   bitops_1.0-7                  rhdf5filters_1.7.0            cachem_1.0.6                  DelayedArray_0.21.2
# [41] assertthat_0.2.1              promises_1.2.0.1              BiocIO_1.5.0                  scales_1.1.1                  beeswarm_0.4.0
# [46] gtable_0.3.0                  beachmat_2.11.0               benchmarkmeData_1.0.4         spam_2.8-0                    rlang_1.0.4
# [51] scatterplot3d_0.3-41          rtracklayer_1.55.4            lazyeval_0.2.2                BiocManager_1.30.16           yaml_2.3.5
# [56] httpuv_1.6.5                  tools_4.2.0                   usethis_2.1.5                 ggplot2_3.3.6                 ellipsis_0.3.2
# [61] jquerylib_0.1.4               RColorBrewer_1.1-3            Rcpp_1.0.8.3                  sparseMatrixStats_1.7.0       progress_1.2.2
# [66] zlibbioc_1.41.0               purrr_0.3.4.9000              RCurl_1.98-1.6                prettyunits_1.1.1             viridis_0.6.2
# [71] cowplot_1.1.1                 ggrepel_0.9.1                 cluster_2.1.3                 fs_1.5.2                      magrittr_2.0.3
# [76] data.table_1.14.2             magick_2.7.3                  pkgload_1.2.4                 hms_1.1.1                     mime_0.12
# [81] xtable_1.8-4                  XML_3.99-0.9                  gridExtra_2.3                 scater_1.23.6                 testthat_3.1.3
# [86] compiler_4.2.0                tibble_3.1.8                  maps_3.4.0                    crayon_1.5.1                  R.oo_1.24.0
# [91] htmltools_0.5.2               later_1.3.0                   tidyr_1.2.0                   DBI_1.1.2                     ExperimentHub_2.3.7
# [96] dbplyr_2.1.1                  rappdirs_0.3.3                Matrix_1.4-1                  brio_1.1.3                    cli_3.2.0
# [101] R.methodsS3_1.8.1             benchmarkme_1.0.7             parallel_4.2.0                metapod_1.3.0                 dotCall64_1.0-1
# [106] igraph_1.3.0                  golem_0.3.2                   pkgconfig_2.0.3               GenomicAlignments_1.31.2      plotly_4.10.0
# [111] xml2_1.3.3                    roxygen2_7.1.2                foreach_1.5.2                 bslib_0.3.1                   vipor_0.4.5
# [116] dqrng_0.3.0                   XVector_0.35.0                attempt_0.3.1                 stringr_1.4.0                 digest_0.6.29
# [121] Biostrings_2.63.3             cellranger_1.1.0              edgeR_3.37.1                  DelayedMatrixStats_1.17.0     restfulr_0.0.13
# [126] curl_4.3.2                    shiny_1.7.1                   Rsamtools_2.11.0              rjson_0.2.21                  lifecycle_1.0.1
# [131] jsonlite_1.8.0                Rhdf5lib_1.17.3               BiocNeighbors_1.13.0          desc_1.4.1                    viridisLite_0.4.0
# [136] limma_3.51.7                  fansi_1.0.3                   pillar_1.7.0                  lattice_0.20-45               KEGGREST_1.35.0
# [141] fastmap_1.1.0                 httr_1.4.2                    interactiveDisplayBase_1.33.0 glue_1.6.2                    png_0.1-7
# [146] iterators_1.0.14              Polychrome_1.3.1              bluster_1.5.1                 BiocVersion_3.15.2            bit_4.0.4
# [151] sass_0.4.1                    stringi_1.7.6                 HDF5Array_1.23.2              blob_1.2.3                    BiocSingular_1.11.0
# [156] AnnotationHub_3.3.11          memoise_2.0.1                 irlba_2.3.5
