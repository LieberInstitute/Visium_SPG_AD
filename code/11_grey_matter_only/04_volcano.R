library("here")
library("sessioninfo")
library("spatialLIBD")
library("EnhancedVolcano")

## Plot output directory
dir_plots <- here::here(
    "plots",
    "11_grey_matter_only",
    "wholegenome"
)
stopifnot(file.exists(dir_plots))

## Locate data directory
dir_rdata <- here::here(
        "code",
        "05_deploy_app_wholegenome"
    )

load(file.path(dir_rdata, "Visium_IF_AD_modeling_results.Rdata"),
    verbose = TRUE
)
sce_pseudo <-
    readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))


## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

subset_sig <- subset(sig_genes, model_type == "enrichment")
tapply(subset_sig$fdr, subset_sig$test, summary)
# $Ab
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00123 0.29067 0.56992 0.54608 0.80964 0.99996
#
# $both
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.7022  0.9996  0.9996  0.9995  0.9996  0.9999
#
# $n_Ab
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0001444 0.0825055 0.2873968 0.3690669 0.6320528 0.9996572
#
# $n_both
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.6274  0.9918  0.9918  0.9924  0.9918  0.9999
#
# $n_pTau
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.9607  0.9607  0.9607  0.9632  0.9607  0.9999
#
# $none
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.972   0.972   0.972   0.974   0.972   1.000
#
# $pTau
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.9327  0.9327  0.9327  0.9378  0.9327  0.9992

tapply(subset_sig$pval, subset_sig$test, summary)
# $Ab
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000001 0.0726954 0.2850369 0.3580220 0.6072802 0.9999650
#
# $both
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0001801 0.3457502 0.5506844 0.5492825 0.7591084 0.9999096
#
# $n_Ab
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.02064 0.14372 0.27219 0.47414 0.99966
#
# $n_both
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000743 0.3677204 0.5505807 0.5568566 0.7560552 0.9999460
#
# $n_pTau
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.03855 0.41855 0.55950 0.57418 0.73158 0.99986
#
# $none
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.002288 0.332950 0.520259 0.529627 0.730310 0.999958
#
# $pTau
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0005831 0.3810386 0.5279638 0.5440502 0.7012146 0.9991950

## Build a data.frame with the info needed first
## Adapted from https://github.com/LieberInstitute/spatial_hpc/blob/62334e88788c2a7010f2d500418639769121ecf1/code/08_pseudobulk/PRECAST/plot_volcano.R#L75-L80

# testname <- "Ab"
make_volcano <- function(testname) {
    subset_sig <- subset(sig_genes, test == testname & model_type == "enrichment")
    df <- data.frame(
        gene = subset_sig$gene,
        logFC = subset_sig$logFC,
        FDR = subset_sig$fdr,
        sig = subset_sig$fdr < 0.1,
        pval = subset_sig$pval
    )

    if(testname == "Ab") {
        selected <- c("UBE2A", "PSMC4", "IDI1", "NINJ1")
    } else if (testname == "n_Ab") {
        selected <- c("C3", "PPP3CA", "UCHL1", "SST")
    } else {
        selected <- df$gene[df$sig]
    }
    n_overlaps <- ifelse(length(selected) > 0, Inf, 15)
    p_cut <- ifelse(min(subset_sig$fdr) < 0.1, subset_sig$pval[subset_sig$fdr < 0.1][which.max(subset_sig$fdr[subset_sig$fdr < 0.1])], 1 / 1e1000)

    ## Adapted from https://github.com/LieberInstitute/spatial_hpc/blob/62334e88788c2a7010f2d500418639769121ecf1/code/08_pseudobulk/PRECAST/plot_volcano.R#L580-L591
    EnhancedVolcano(df,
        lab = df$gene,
        x = 'logFC',
        y = 'pval',
        FCcutoff = 0,
        pCutoff = p_cut,
        ylab = "-log10 pval",
        legendLabels = c('Not sig.','Not sig.', 'FDR < 0.1',
            'FDR < 0.1'),
        col = c("grey30", "grey30", "royalblue", "red2"),
        title = paste(testname, "> rest"),
        subtitle = "",
        selectLab = selected,
        max.overlaps = n_overlaps,
        drawConnectors = TRUE,
        lengthConnectors = unit(0.013, "npc"),
        directionConnectors = "y",
        caption = paste0("total = ", nrow(df), " genes")
    )
}

pdf(file.path(dir_plots, "volcano_plots.pdf"), width = 5)
for(i in sort(unique(subset_sig$test))) {
    message(Sys.time(), " making volcano plot for ", i)
    print(make_volcano(i))
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 (2022-10-31)
#  os       macOS Ventura 13.0.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/Mexico_City
#  date     2023-02-24
#  rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  AnnotationDbi            1.60.0    2022-11-01 [1] Bioconductor
#  AnnotationHub            3.6.0     2022-11-01 [1] Bioconductor
#  assertthat               0.2.1     2019-03-21 [1] CRAN (R 4.2.0)
#  attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.2.0)
#  beachmat                 2.14.0    2022-11-01 [1] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.2.0)
#  benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.2.0)
#  benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.2.0)
#  Biobase                * 2.58.0    2022-11-01 [1] Bioconductor
#  BiocFileCache            2.6.0     2022-11-01 [1] Bioconductor
#  BiocGenerics           * 0.44.0    2022-11-01 [1] Bioconductor
#  BiocIO                   1.8.0     2022-11-01 [1] Bioconductor
#  BiocManager              1.30.19   2022-10-25 [1] CRAN (R 4.2.0)
#  BiocNeighbors            1.16.0    2022-11-01 [1] Bioconductor
#  BiocParallel             1.32.5    2022-12-25 [1] Bioconductor
#  BiocSingular             1.14.0    2022-11-01 [1] Bioconductor
#  BiocVersion              3.16.0    2022-09-20 [1] Bioconductor
#  Biostrings               2.66.0    2022-11-01 [1] Bioconductor
#  bit                      4.0.5     2022-11-15 [1] CRAN (R 4.2.2)
#  bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.2.0)
#  bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.2.0)
#  blob                     1.2.3     2022-04-10 [1] CRAN (R 4.2.0)
#  brio                     1.1.3     2021-11-30 [1] CRAN (R 4.2.0)
#  bslib                    0.4.2     2022-12-16 [1] CRAN (R 4.2.2)
#  cachem                   1.0.6     2021-08-19 [1] CRAN (R 4.2.0)
#  callr                    3.7.3     2022-11-02 [1] CRAN (R 4.2.2)
#  cli                      3.6.0     2023-01-09 [1] CRAN (R 4.2.0)
#  codetools                0.2-19    2023-02-01 [1] CRAN (R 4.2.0)
#  colorout                 1.2-2     2022-03-01 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.2.0)
#  config                   0.3.1     2020-12-17 [1] CRAN (R 4.2.0)
#  cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.2.0)
#  crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.2.0)
#  curl                     5.0.0     2023-01-12 [1] CRAN (R 4.2.0)
#  data.table               1.14.8    2023-02-17 [1] CRAN (R 4.2.2)
#  DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.2.0)
#  dbplyr                   2.3.0     2023-01-16 [1] CRAN (R 4.2.0)
#  DelayedArray             0.24.0    2022-11-01 [1] Bioconductor
#  DelayedMatrixStats       1.20.0    2022-11-01 [1] Bioconductor
#  desc                     1.4.2     2022-09-08 [1] CRAN (R 4.2.0)
#  devtools               * 2.4.5     2022-10-11 [1] CRAN (R 4.2.0)
#  digest                   0.6.31    2022-12-11 [1] CRAN (R 4.2.0)
#  doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.2.0)
#  dotCall64                1.0-2     2022-10-03 [1] CRAN (R 4.2.1)
#  dplyr                    1.1.0     2023-01-29 [1] CRAN (R 4.2.0)
#  dqrng                    0.3.0     2021-05-01 [1] CRAN (R 4.2.0)
#  DropletUtils             1.18.1    2022-11-23 [1] Bioconductor
#  DT                       0.27      2023-01-17 [1] CRAN (R 4.2.0)
#  edgeR                    3.40.2    2023-01-22 [1] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.2.0)
#  EnhancedVolcano        * 1.16.0    2022-11-07 [1] Bioconductor
#  ExperimentHub            2.6.0     2022-11-01 [1] Bioconductor
#  fansi                    1.0.4     2023-01-22 [1] CRAN (R 4.2.0)
#  farver                   2.1.1     2022-07-06 [1] CRAN (R 4.2.1)
#  fastmap                  1.1.0     2021-01-25 [1] CRAN (R 4.2.0)
#  fields                   14.1      2022-08-12 [1] CRAN (R 4.2.0)
#  filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.2.0)
#  foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.2.0)
#  fs                       1.6.1     2023-02-06 [1] CRAN (R 4.2.0)
#  generics                 0.1.3     2022-07-05 [1] CRAN (R 4.2.0)
#  GenomeInfoDb           * 1.34.9    2023-02-02 [1] Bioconductor
#  GenomeInfoDbData         1.2.9     2022-11-02 [1] Bioconductor
#  GenomicAlignments        1.34.0    2022-11-01 [1] Bioconductor
#  GenomicRanges          * 1.50.2    2022-12-18 [1] Bioconductor
#  ggbeeswarm               0.7.1     2022-12-16 [1] CRAN (R 4.2.2)
#  ggplot2                * 3.4.1     2023-02-10 [1] CRAN (R 4.2.0)
#  ggrepel                * 0.9.3     2023-02-03 [1] CRAN (R 4.2.0)
#  glue                     1.6.2     2022-02-24 [1] CRAN (R 4.2.0)
#  golem                    0.3.5     2022-10-18 [1] CRAN (R 4.2.0)
#  gridExtra                2.3       2017-09-09 [1] CRAN (R 4.2.0)
#  gtable                   0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
#  HDF5Array                1.26.0    2022-11-01 [1] Bioconductor
#  here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.2.0)
#  hms                      1.1.2     2022-08-19 [1] CRAN (R 4.2.0)
#  htmltools                0.5.4     2022-12-07 [1] CRAN (R 4.2.0)
#  htmlwidgets              1.6.1     2023-01-07 [1] CRAN (R 4.2.0)
#  httpuv                   1.6.9     2023-02-14 [1] CRAN (R 4.2.0)
#  httr                     1.4.4     2022-08-17 [1] CRAN (R 4.2.0)
#  interactiveDisplayBase   1.36.0    2022-11-01 [1] Bioconductor
#  IRanges                * 2.32.0    2022-11-01 [1] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [1] CRAN (R 4.2.0)
#  jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.2.0)
#  jsonlite                 1.8.4     2022-12-06 [1] CRAN (R 4.2.0)
#  KEGGREST                 1.38.0    2022-11-01 [1] Bioconductor
#  knitr                    1.42      2023-01-25 [1] CRAN (R 4.2.0)
#  labeling                 0.4.2     2020-10-20 [1] CRAN (R 4.2.0)
#  later                    1.3.0     2021-08-18 [1] CRAN (R 4.2.0)
#  lattice                  0.20-45   2021-09-22 [1] CRAN (R 4.2.2)
#  lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.2.0)
#  lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
#  limma                    3.54.1    2023-01-26 [1] Bioconductor
#  locfit                   1.5-9.7   2023-01-02 [1] CRAN (R 4.2.0)
#  lubridate                1.9.2     2023-02-10 [1] CRAN (R 4.2.0)
#  magick                   2.7.3     2021-08-18 [1] CRAN (R 4.2.0)
#  magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.2.0)
#  maps                     3.4.1     2022-10-30 [1] CRAN (R 4.2.0)
#  Matrix                   1.5-3     2022-11-11 [1] CRAN (R 4.2.0)
#  MatrixGenerics         * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats            * 0.63.0    2022-11-18 [1] CRAN (R 4.2.0)
#  memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.2.0)
#  mime                     0.12      2021-09-28 [1] CRAN (R 4.2.0)
#  miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.2.0)
#  munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.2.0)
#  paletteer                1.5.0     2022-10-19 [1] CRAN (R 4.2.0)
#  pillar                   1.8.1     2022-08-19 [1] CRAN (R 4.2.0)
#  pkgbuild                 1.4.0     2022-11-27 [1] CRAN (R 4.2.2)
#  pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.2.0)
#  pkgload                  1.3.2     2022-11-16 [1] CRAN (R 4.2.2)
#  plotly                   4.10.1    2022-11-07 [1] CRAN (R 4.2.0)
#  png                      0.1-8     2022-11-29 [1] CRAN (R 4.2.0)
#  prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.2.0)
#  processx                 3.8.0     2022-10-26 [1] CRAN (R 4.2.0)
#  profvis                  0.3.7     2020-11-02 [1] CRAN (R 4.2.0)
#  promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.2.0)
#  prompt                   1.0.1     2022-03-01 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                       1.7.2     2022-10-26 [1] CRAN (R 4.2.0)
#  purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.2.0)
#  R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
#  R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
#  R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.2.0)
#  R6                       2.5.1     2021-08-19 [1] CRAN (R 4.2.0)
#  rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.2.0)
#  RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.2.0)
#  Rcpp                     1.0.10    2023-01-22 [1] CRAN (R 4.2.0)
#  RCurl                    1.98-1.10 2023-01-27 [1] CRAN (R 4.2.0)
#  rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.2.0)
#  remotes                  2.4.2     2021-11-30 [1] CRAN (R 4.2.0)
#  restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.2.0)
#  rhdf5                    2.42.0    2022-11-01 [1] Bioconductor
#  rhdf5filters             1.10.0    2022-11-01 [1] Bioconductor
#  Rhdf5lib                 1.20.0    2022-11-01 [1] Bioconductor
#  rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.2.0)
#  rlang                    1.0.6     2022-09-24 [1] CRAN (R 4.2.0)
#  roxygen2                 7.2.3     2022-12-08 [1] CRAN (R 4.2.0)
#  rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.2.0)
#  Rsamtools                2.14.0    2022-11-01 [1] Bioconductor
#  RSQLite                  2.3.0     2023-02-17 [1] CRAN (R 4.2.2)
#  rsthemes                 0.3.1     2022-03-01 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi               0.14      2022-08-22 [1] CRAN (R 4.2.0)
#  rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.2.0)
#  rtracklayer              1.58.0    2022-11-01 [1] Bioconductor
#  S4Vectors              * 0.36.1    2022-12-07 [1] Bioconductor
#  sass                     0.4.5     2023-01-24 [1] CRAN (R 4.2.0)
#  ScaledMatrix             1.6.0     2022-11-01 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [1] CRAN (R 4.2.0)
#  scater                   1.26.1    2022-11-13 [1] Bioconductor
#  scuttle                  1.8.4     2023-01-22 [1] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.2.0)
#  shiny                    1.7.4     2022-12-15 [1] CRAN (R 4.2.2)
#  shinyWidgets             0.7.6     2023-01-08 [1] CRAN (R 4.2.0)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [1] CRAN (R 4.2.0)
#  sparseMatrixStats        1.10.0    2022-11-01 [1] Bioconductor
#  SpatialExperiment      * 1.8.0     2022-11-01 [1] Bioconductor
#  spatialLIBD            * 1.11.7    2023-02-17 [1] Github (LieberInstitute/spatialLIBD@bce82bf)
#  statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.2.0)
#  stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.2.0)
#  stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.2.0)
#  SummarizedExperiment   * 1.28.0    2022-11-01 [1] Bioconductor
#  suncalc                  0.5.1     2022-09-29 [1] CRAN (R 4.2.0)
#  testthat               * 3.1.6     2022-12-09 [1] CRAN (R 4.2.0)
#  tibble                   3.1.8     2022-07-22 [1] CRAN (R 4.2.1)
#  tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.2.0)
#  tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.2.0)
#  timechange               0.2.0     2023-01-11 [1] CRAN (R 4.2.0)
#  urlchecker               1.0.1     2021-11-30 [1] CRAN (R 4.2.0)
#  usethis                * 2.1.6     2022-05-25 [1] CRAN (R 4.2.0)
#  utf8                     1.2.3     2023-01-31 [1] CRAN (R 4.2.0)
#  vctrs                    0.5.2     2023-01-23 [1] CRAN (R 4.2.0)
#  vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.2.0)
#  viridis                  0.6.2     2021-10-13 [1] CRAN (R 4.2.0)
#  viridisLite              0.4.1     2022-08-22 [1] CRAN (R 4.2.0)
#  withr                    2.5.0     2022-03-03 [1] CRAN (R 4.2.0)
#  xfun                     0.37      2023-01-31 [1] CRAN (R 4.2.0)
#  XML                      3.99-0.13 2022-12-04 [1] CRAN (R 4.2.0)
#  xml2                     1.3.3     2021-11-30 [1] CRAN (R 4.2.0)
#  xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.2.0)
#  XVector                  0.38.0    2022-11-01 [1] Bioconductor
#  yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.2.0)
#  zlibbioc                 1.44.0    2022-11-01 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
