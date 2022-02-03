library("spatialLIBD")
library("Polychrome")
library("here")
library("sessioninfo")

plot_SNN_k10 <- function(suffix) {
    spe <- readRDS(
        here::here(
            "processed-data", "08_harmony_BayesSpace", suffix,
            paste0("spe_harmony_", suffix, ".rds")
        )
    )

    cols <- Polychrome::palette36.colors(length(unique(spe$SNN_k10)))
    names(cols) <- sort(unique(spe$SNN_k10))

    vis_grid_clus(
        spe = spe,
        clustervar = "SNN_k10",
        pdf_file = here::here("plots", "08_harmony_BayesSpace", suffix, "graphbased_SNN_k10.pdf"),
        sort_clust = FALSE,
        colors = cols,
        spatial = FALSE,
        point_size = 2
    )

    pdf(here::here("plots", "08_harmony_BayesSpace", suffix, "graphbased_SNN_k10_cut_at.pdf"), height = 24, width = 36)
    for (k in 4:28) {
        k_nice <- sprintf("%02d", k)

        cols <- Polychrome::palette36.colors(k)
        k_var <- paste0("SNN_k10_k", k)
        names(cols) <- sort(unique(colData(spe)[[k_var]]))

        p_list <- vis_grid_clus(
            spe = spe,
            clustervar = k_var,
            return_plots = TRUE,
            sort_clust = FALSE,
            colors = cols,
            spatial = FALSE,
            point_size = 2
        )
        print(cowplot::plot_grid(plotlist = p_list))
    }
    dev.off()

    return(NULL)
}

plot_SNN_k10("wholegenome")
plot_SNN_k10("targeted")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.2 Patched (2021-11-04 r81138)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-02-03
#  pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date (UTC) lib source
#  AnnotationDbi            1.56.2   2021-11-09 [2] Bioconductor
#  AnnotationHub            3.2.1    2022-01-23 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.2)
#  beachmat                 2.10.0   2021-10-26 [2] Bioconductor
#  beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.2)
#  benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.2)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.2)
#  Biobase                * 2.54.0   2021-10-26 [2] Bioconductor
#  BiocFileCache            2.2.1    2022-01-23 [2] Bioconductor
#  BiocGenerics           * 0.40.0   2021-10-26 [2] Bioconductor
#  BiocIO                   1.4.0    2021-10-26 [2] Bioconductor
#  BiocManager              1.30.16  2021-06-15 [2] CRAN (R 4.1.2)
#  BiocNeighbors            1.12.0   2021-10-26 [1] Bioconductor
#  BiocParallel             1.28.3   2021-12-09 [2] Bioconductor
#  BiocSingular             1.10.0   2021-10-26 [1] Bioconductor
#  BiocVersion              3.14.0   2021-05-19 [2] Bioconductor
#  Biostrings               2.62.0   2021-10-26 [2] Bioconductor
#  bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  blob                     1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
#  brio                     1.1.3    2021-11-30 [2] CRAN (R 4.1.2)
#  bslib                    0.3.1    2021-10-06 [2] CRAN (R 4.1.2)
#  cachem                   1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
#  callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
#  cli                      3.1.1    2022-01-20 [2] CRAN (R 4.1.2)
#  codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
#  colorout                 1.2-2    2021-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.2)
#  cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.2)
#  crayon                   1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
#  curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
#  data.table               1.14.2   2021-09-27 [2] CRAN (R 4.1.2)
#  DBI                      1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
#  dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
#  DelayedArray             0.20.0   2021-10-26 [2] Bioconductor
#  DelayedMatrixStats       1.16.0   2021-10-26 [2] Bioconductor
#  desc                     1.4.0    2021-09-28 [2] CRAN (R 4.1.2)
#  digest                   0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
#  dockerfiler              0.1.4    2021-09-03 [1] CRAN (R 4.1.2)
#  doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)
#  dplyr                    1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.2)
#  DropletUtils             1.14.2   2022-01-09 [1] Bioconductor
#  DT                       0.20     2021-11-15 [2] CRAN (R 4.1.2)
#  edgeR                    3.36.0   2021-10-26 [2] Bioconductor
#  ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  ExperimentHub            2.2.1    2022-01-23 [2] Bioconductor
#  fansi                    1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
#  farver                   2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   13.3     2021-10-30 [2] CRAN (R 4.1.2)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.2    2022-02-02 [2] CRAN (R 4.1.2)
#  fs                       1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
#  generics                 0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
#  GenomeInfoDb           * 1.30.1   2022-01-30 [2] Bioconductor
#  GenomeInfoDbData         1.2.7    2021-11-01 [2] Bioconductor
#  GenomicAlignments        1.30.0   2021-10-26 [2] Bioconductor
#  GenomicRanges          * 1.46.1   2021-11-18 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.2)
#  ggplot2                  3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  glue                     1.6.1    2022-01-22 [2] CRAN (R 4.1.2)
#  golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.2)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.22.1   2021-11-14 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.2)
#  htmltools                0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
#  htmlwidgets              1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
#  httpuv                   1.6.5    2022-01-05 [2] CRAN (R 4.1.2)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  interactiveDisplayBase   1.32.0   2021-10-26 [2] Bioconductor
#  IRanges                * 2.28.0   2021-10-26 [2] Bioconductor
#  irlba                    2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.3    2022-01-17 [2] CRAN (R 4.1.2)
#  KEGGREST                 1.34.0   2021-10-26 [2] Bioconductor
#  knitr                    1.37     2021-12-16 [2] CRAN (R 4.1.2)
#  labeling                 0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
#  later                    1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
#  lattice                  0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
#  limma                    3.50.0   2021-10-26 [2] Bioconductor
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
#  magrittr                 2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
#  maps                     3.4.0    2021-09-25 [2] CRAN (R 4.1.2)
#  Matrix                   1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
#  MatrixGenerics         * 1.6.0    2021-10-26 [2] Bioconductor
#  matrixStats            * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
#  memoise                  2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
#  mime                     0.12     2021-09-28 [2] CRAN (R 4.1.2)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
#  pkgbuild                 1.3.1    2021-12-20 [2] CRAN (R 4.1.2)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.2.4    2021-11-30 [2] CRAN (R 4.1.2)
#  plotly                   4.10.0   2021-10-09 [2] CRAN (R 4.1.2)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome             * 1.3.1    2021-07-16 [1] CRAN (R 4.1.2)
#  prettyunits              1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
#  processx                 3.5.2    2021-04-30 [2] CRAN (R 4.1.0)
#  promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  ps                       1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
#  purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                  2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
#  R6                       2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
#  rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                     1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
#  RCurl                    1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
#  remotes                  2.4.2    2021-11-30 [2] CRAN (R 4.1.2)
#  restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                    2.38.0   2021-10-26 [2] Bioconductor
#  rhdf5filters             1.6.0    2021-10-26 [2] Bioconductor
#  Rhdf5lib                 1.16.0   2021-10-26 [2] Bioconductor
#  rjson                    0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
#  rlang                    1.0.0    2022-01-26 [2] CRAN (R 4.1.2)
#  rmote                    0.3.4    2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.1.2    2021-09-08 [2] CRAN (R 4.1.2)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools                2.10.0   2021-10-26 [2] Bioconductor
#  RSQLite                  2.2.9    2021-12-06 [2] CRAN (R 4.1.2)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.2)
#  rtracklayer              1.54.0   2021-10-26 [2] Bioconductor
#  S4Vectors              * 0.32.3   2021-11-21 [2] Bioconductor
#  sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
#  ScaledMatrix             1.2.0    2021-10-26 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                   1.22.0   2021-10-26 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.2)
#  scuttle                  1.4.0    2021-10-26 [1] Bioconductor
#  servr                    0.24     2021-11-16 [1] CRAN (R 4.1.2)
#  sessioninfo            * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
#  shiny                    1.7.1    2021-10-02 [2] CRAN (R 4.1.2)
#  shinyWidgets             0.6.3    2022-01-10 [1] CRAN (R 4.1.2)
#  SingleCellExperiment   * 1.16.0   2021-10-26 [2] Bioconductor
#  spam                     2.8-0    2022-01-06 [2] CRAN (R 4.1.2)
#  sparseMatrixStats        1.6.0    2021-10-26 [2] Bioconductor
#  SpatialExperiment      * 1.4.0    2021-10-26 [1] Bioconductor
#  spatialLIBD            * 1.6.5    2022-01-12 [1] Bioconductor
#  stringi                  1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.24.0   2021-10-26 [2] Bioconductor
#  testthat                 3.1.2    2022-01-20 [2] CRAN (R 4.1.2)
#  tibble                   3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
#  tidyr                    1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
#  tidyselect               1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  usethis                  2.1.5    2021-12-09 [2] CRAN (R 4.1.2)
#  utf8                     1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.2)
#  viridis                  0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
#  viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  withr                    2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
#  xfun                     0.29     2021-12-14 [2] CRAN (R 4.1.2)
#  XML                      3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
#  xml2                     1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
#  xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                  0.34.0   2021-10-26 [2] Bioconductor
#  yaml                     2.2.2    2022-01-25 [2] CRAN (R 4.1.2)
#  zlibbioc                 1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────
