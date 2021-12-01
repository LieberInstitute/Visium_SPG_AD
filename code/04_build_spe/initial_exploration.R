library("SpatialExperiment")
library("here")
library("spatialLIBD")
library("readxl")
library("RColorBrewer")
library("sessioninfo")

load(here("processed-data", "spe", "spe.Rdata"), verbose = TRUE)
load(here("processed-data", "spe", "spe_targeted.Rdata"), verbose = TRUE)

dir.create(here("plots", "initial_exploration"), showWarnings = FALSE)

## Define order of samples for the grid plots
slide_order <- c("V10A27106", "V10T31036", "V10A27004")
sample_order <- unlist(sapply(slide_order, function(i) {
    sort(unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
sample_order
#            V10A271061            V10A271062            V10A271063
# "V10A27106_A1_Br3874" "V10A27106_B1_Br3854" "V10A27106_C1_Br3873"
#            V10A271064            V10T310361            V10T310362
# "V10A27106_D1_Br3880" "V10T31036_A1_Br3874" "V10T31036_B1_Br3854"
#            V10T310363            V10T310364            V10A270041
# "V10T31036_C1_Br3873" "V10T31036_D1_Br3880" "V10A27004_A1_Br3874"
#            V10A270042
# "V10A27004_D1_Br3880"

## Check the segmentation information
segmentation_variables <-
    c(
        "NAbeta",
        "PAbeta",
        "NDAPI",
        "PDAPI",
        "NGFAP",
        "PGFAP",
        "NLipofuscin",
        "PLipofuscin",
        "NMAP2",
        "PMAP2",
        "NpTau",
        "PpTau"
    )

for(seg_var in segmentation_variables) {

    seg_grid <-
        vis_grid_gene(
            spe,
            geneid = seg_var,
            return_plots = TRUE,
            spatial = FALSE,
            cont_colors = viridisLite::magma(21, direction = -1),
            minCount = -1
        )
    pdf(
        here(
            "plots",
            "initial_exploration",
            paste0("segmentation_info_", seg_var, ".pdf")
        ),
        height = 24,
        width = 36
    )
    print(cowplot::plot_grid(plotlist = seg_grid[sample_order]))
    dev.off()
}

#### the code below is no longer needed since we now have segmentation data for
#### all images
## Check the segmantation preliminary results on one sample
# pdf(here("plots", "initial_exploration", "segmentation_info.pdf"), height = 8, width = 9)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NpTau", spatial = FALSE, cont_colors = viridisLite::mako(21, direction = -1))
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PpTau", spatial = FALSE, cont_colors = viridisLite::mako(21, direction = -1))
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NAbeta", spatial = FALSE, cont_colors = viridisLite::mako(21, direction = -1))
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PAbeta", spatial = FALSE, cont_colors = viridisLite::mako(21, direction = -1))
#
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NpTau", spatial = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PpTau", spatial = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NAbeta", spatial = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PAbeta", spatial = FALSE)
#
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NpTau", spatial = TRUE, viridis = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PpTau", spatial = TRUE, viridis = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "NAbeta", spatial = TRUE, viridis = FALSE)
# vis_gene(spe, sampleid = "V10A27106_D1_Br3880", geneid = "PAbeta", spatial = TRUE, viridis = FALSE)
# dev.off()


## Make grids for the GraphBased cluster results
length(unique(spe$`10x_graphclust`))
# [1] 9
cols <- RColorBrewer::brewer.pal(length(unique(spe$`10x_graphclust`)), "Set1")
names(cols) <- seq_len(length(cols))
clus_list <- vis_grid_clus(spe,
    clustervar = "10x_graphclust",
    pdf_file = NULL,
    sort_clust = TRUE,
    colors = cols,
    return_plots = TRUE,
    spatial = FALSE
)

pdf(here("plots", "initial_exploration", "wholegenome_graph_based.pdf"), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = clus_list[sample_order]))
dev.off()

length(unique(spe_targeted$`10x_graphclust`))
# [1] 6
cols_targeted <- RColorBrewer::brewer.pal(length(unique(spe_targeted$`10x_graphclust`)), "Dark2")
names(cols_targeted) <- seq_len(length(cols_targeted))
clus_list_targeted <- vis_grid_clus(spe_targeted,
    clustervar = "10x_graphclust",
    pdf_file = NULL,
    sort_clust = TRUE,
    colors = cols_targeted,
    return_plots = TRUE,
    spatial = FALSE
)
pdf(here("plots", "initial_exploration", "targeted_graph_based.pdf"), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = clus_list_targeted[sample_order]))
dev.off()


## Read in list of AD genes
ad_genes_raw <- read_xlsx(here("raw-data", "10X_NS targeted gene_AD genes.xlsx"))
ad_genes <- paste0(ad_genes_raw[[2]], "; ", ad_genes_raw[[1]])

addmargins(table(
    "WholeGenome" = ad_genes %in% rowRanges(spe)$gene_search,
    "TargetedSequencing" = ad_genes %in% rowRanges(spe_targeted)$gene_search
))
#            TargetedSequencing
# WholeGenome FALSE TRUE Sum
#       FALSE     3    2   5
#       TRUE      0  134 134
#       Sum       3  136 139
ad_genes[!ad_genes %in% rowRanges(spe)$gene_search]
# [1] "ENSG00000070748; CHAT"    "ENSG00000111537; IFNG"    "ENSG00000113520; IL4"     "ENSG00000254647; INS"
# [5] "ENSG00000187714; SLC18A3"
ad_genes[!ad_genes %in% rowRanges(spe_targeted)$gene_search]
# [1] "ENSG00000070748; CHAT" "ENSG00000111537; IFNG" "ENSG00000113520; IL4"

make_gene_grids <- function(spatial) {
    genes <- ad_genes[ad_genes %in% rowRanges(spatial)$gene_search]
    for(g in genes) {
        print(vis_gene(
            spatial,
            sampleid = "V10A27106_D1_Br3880",
            geneid = g,
            spatial = FALSE,
            assayname = "counts",
            viridis = FALSE
        ))
    }
    return(NULL)
}

pdf(here("plots", "initial_exploration", "wholegenome_V10A27106_D1_Br3880_AD_genes.pdf"), height = 8, width = 9)
make_gene_grids(spe)
dev.off()

pdf(here("plots", "initial_exploration", "targeted_V10A27106_D1_Br3880_AD_genes.pdf"), height = 8, width = 9)
make_gene_grids(spe_targeted)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.1 Patched (2021-08-13 r80752)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-10-12
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date       lib source
#  AnnotationDbi            1.55.1   2021-06-07 [2] Bioconductor
#  AnnotationHub            3.1.5    2021-08-12 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.1)
#  beachmat                 2.9.1    2021-08-11 [2] Bioconductor
#  beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.1)
#  benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.1)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.1)
#  Biobase                * 2.53.0   2021-05-19 [2] Bioconductor
#  BiocFileCache            2.1.1    2021-06-23 [2] Bioconductor
#  BiocGenerics           * 0.39.2   2021-08-18 [1] Bioconductor
#  BiocManager              1.30.16  2021-06-15 [2] CRAN (R 4.1.0)
#  BiocNeighbors            1.11.0   2021-05-19 [1] Bioconductor
#  BiocParallel             1.27.17  2021-10-10 [1] Bioconductor
#  BiocSingular             1.9.1    2021-06-08 [1] Bioconductor
#  BiocVersion              3.14.0   2021-05-19 [2] Bioconductor
#  Biostrings               2.61.2   2021-08-04 [2] Bioconductor
#  bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  blob                     1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
#  bslib                    0.2.5.1  2021-05-18 [2] CRAN (R 4.1.0)
#  cachem                   1.0.6    2021-08-19 [1] CRAN (R 4.1.1)
#  callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
#  cellranger               1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
#  cli                      3.0.1    2021-07-17 [2] CRAN (R 4.1.0)
#  codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.1)
#  colorout                 1.2-2    2021-09-13 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.1)
#  cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.1)
#  crayon                   1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
#  data.table               1.14.0   2021-02-21 [2] CRAN (R 4.1.0)
#  DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
#  DelayedArray             0.19.4   2021-09-23 [1] Bioconductor
#  DelayedMatrixStats       1.15.2   2021-08-05 [2] Bioconductor
#  desc                     1.3.0    2021-03-05 [2] CRAN (R 4.1.0)
#  digest                   0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  dockerfiler              0.1.4    2021-09-03 [1] CRAN (R 4.1.1)
#  doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)
#  dplyr                    1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.1)
#  DropletUtils             1.13.4   2021-09-19 [1] Bioconductor
#  DT                       0.18     2021-04-14 [2] CRAN (R 4.1.0)
#  edgeR                    3.35.0   2021-05-19 [2] Bioconductor
#  ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  ExperimentHub            2.1.4    2021-07-27 [2] Bioconductor
#  fansi                    0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   12.5     2021-06-25 [2] CRAN (R 4.1.0)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics                 0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb           * 1.29.8   2021-09-05 [1] Bioconductor
#  GenomeInfoDbData         1.2.6    2021-05-21 [2] Bioconductor
#  GenomicRanges          * 1.45.0   2021-05-19 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.1)
#  ggplot2                  3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  glue                     1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.1)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.21.0   2021-05-19 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.1)
#  htmltools                0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets              1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                   1.6.1    2021-05-07 [2] CRAN (R 4.1.0)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  interactiveDisplayBase   1.31.2   2021-07-30 [2] Bioconductor
#  IRanges                * 2.27.2   2021-08-18 [1] Bioconductor
#  irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  KEGGREST                 1.33.0   2021-05-19 [2] Bioconductor
#  knitr                    1.33     2021-04-24 [2] CRAN (R 4.1.0)
#  later                    1.2.0    2021-04-23 [2] CRAN (R 4.1.0)
#  lattice                  0.20-44  2021-05-02 [3] CRAN (R 4.1.1)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  limma                    3.49.4   2021-08-08 [2] Bioconductor
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.7.2    2021-05-02 [2] CRAN (R 4.1.0)
#  magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  maps                     3.3.0    2018-04-03 [2] CRAN (R 4.1.0)
#  Matrix                   1.3-4    2021-06-01 [3] CRAN (R 4.1.1)
#  MatrixGenerics         * 1.5.4    2021-08-26 [1] Bioconductor
#  matrixStats            * 0.61.0   2021-09-17 [1] CRAN (R 4.1.1)
#  memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
#  mime                     0.11     2021-06-23 [2] CRAN (R 4.1.0)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.6.2    2021-07-29 [2] CRAN (R 4.1.0)
#  pkgbuild                 1.2.0    2020-12-15 [2] CRAN (R 4.1.0)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.2.1    2021-04-06 [2] CRAN (R 4.1.0)
#  plotly                   4.9.4.1  2021-06-18 [2] CRAN (R 4.1.0)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome               1.3.1    2021-07-16 [1] CRAN (R 4.1.1)
#  prettyunits              1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
#  processx                 3.5.2    2021-04-30 [2] CRAN (R 4.1.0)
#  promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  ps                       1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
#  purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                  2.10.1   2020-08-26 [2] CRAN (R 4.1.0)
#  R6                       2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
#  rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  RColorBrewer           * 1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                     1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                    1.98-1.5 2021-09-17 [1] CRAN (R 4.1.1)
#  readxl                 * 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
#  remotes                  2.4.0    2021-06-02 [2] CRAN (R 4.1.0)
#  rhdf5                    2.37.0   2021-05-19 [2] Bioconductor
#  rhdf5filters             1.5.0    2021-05-19 [2] Bioconductor
#  Rhdf5lib                 1.15.2   2021-07-01 [2] Bioconductor
#  rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                    0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
#  rmote                    0.3.4    2021-09-13 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.1.1    2020-06-27 [2] CRAN (R 4.1.0)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  RSQLite                  2.2.8    2021-08-21 [1] CRAN (R 4.1.1)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.1)
#  S4Vectors              * 0.31.5   2021-10-01 [1] Bioconductor
#  sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
#  ScaledMatrix             1.1.0    2021-05-19 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                   1.21.8   2021-10-08 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.1)
#  scuttle                  1.3.1    2021-08-05 [1] Bioconductor
#  servr                    0.23     2021-08-11 [1] CRAN (R 4.1.1)
#  sessioninfo            * 1.1.1    2018-11-05 [1] CRAN (R 4.1.1)
#  shiny                    1.6.0    2021-01-25 [2] CRAN (R 4.1.0)
#  shinyWidgets             0.6.2    2021-09-17 [1] CRAN (R 4.1.1)
#  SingleCellExperiment   * 1.15.1   2021-05-21 [2] Bioconductor
#  spam                     2.7-0    2021-06-25 [2] CRAN (R 4.1.0)
#  sparseMatrixStats        1.5.2    2021-08-05 [2] Bioconductor
#  SpatialExperiment      * 1.3.4    2021-08-24 [1] Bioconductor
#  spatialLIBD            * 1.5.3    2021-08-07 [1] Bioconductor
#  stringi                  1.7.3    2021-07-16 [2] CRAN (R 4.1.0)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.23.5   2021-10-05 [1] Bioconductor
#  testthat                 3.0.4    2021-07-01 [2] CRAN (R 4.1.0)
#  tibble                   3.1.4    2021-08-25 [1] CRAN (R 4.1.1)
#  tidyr                    1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
#  tidyselect               1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  usethis                  2.0.1    2021-02-10 [2] CRAN (R 4.1.0)
#  utf8                     1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.1)
#  viridis                  0.6.1    2021-05-11 [2] CRAN (R 4.1.0)
#  viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  withr                    2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  xfun                     0.25     2021-08-06 [2] CRAN (R 4.1.1)
#  xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
#  xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                  0.33.0   2021-05-19 [2] Bioconductor
#  yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc                 1.39.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lcollado/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
