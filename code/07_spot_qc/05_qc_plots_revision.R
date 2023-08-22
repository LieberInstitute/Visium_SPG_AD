# sgejobs::job_single(
#     "qc_metrics_and_segmentation",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "80G",
#     command = "Rscript 01_qc_metrics_and_segmentation.R",
#     create_logdir = TRUE
# )

library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")

## Load basic SPE data
spe_wholegenome <- readRDS(
    here::here(
        "processed-data", "04_build_spe", "spe_wholegenome.rds"
    )
)
spe_targeted <- readRDS(
    here::here(
        "processed-data", "04_build_spe", "spe_targeted.rds"
    )
)

## Create output directories
dir_plots <- here::here("plots", "07_spot_qc")
dir_rdata <- here::here("processed-data", "07_spot_qc")
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

slide_order <- c("V10A27106", "V10T31036", "V10A27004")
sample_order <- unlist(sapply(slide_order, function(i) {
    sort(unique(spe_wholegenome$sample_id)[grepl(i, unique(spe_wholegenome$sample_id))])
}))
sample_order

## Re-order samples
new_order <- unlist(lapply(sample_order, function(i) {
    which(spe_wholegenome$sample_id == i)
}))
stopifnot(all(seq_len(ncol(spe_wholegenome)) %in% new_order))
spe_wholegenome <- spe_wholegenome[, new_order]
spe_targeted <- spe_targeted[, new_order]

shorten_names <- function(x) {
    gsub(slide_order[3], "S3", gsub(slide_order[2], "S2", gsub(slide_order[1], "S1", x)))
}

spe_wholegenome$sample_id_short <- factor(shorten_names(spe_wholegenome$sample_id), levels = shorten_names(sample_order))
spe_targeted$sample_id_short <- spe_wholegenome$sample_id_short

palette_colors <- paste0(
    rep(c("darkorchid", "sienna", "steelblue"), c(4, 4, 2)),
    c(1:4, 1:4, 3:4)
)
names(palette_colors) <- levels(spe_wholegenome$sample_id_short)


## Re-order samples to match order in other figures
sample_order <- c(
    "V10A27106_B1_Br3854",
    "V10A27106_C1_Br3873",
    "V10A27106_D1_Br3880",
    "V10T31036_B1_Br3854",
    "V10T31036_C1_Br3873",
    "V10T31036_D1_Br3880",
    "V10A27004_D1_Br3880",
    "V10A27106_A1_Br3874",
    "V10T31036_A1_Br3874",
    "V10A27004_A1_Br3874"
)
stopifnot(length(unique(sample_order)) == 10)
new_order <- unlist(lapply(sample_order, function(i) {
    which(spe_wholegenome$sample_id == i)
}))
stopifnot(all(seq_len(ncol(spe_wholegenome)) %in% new_order))
spe_wholegenome <- spe_wholegenome[, new_order]
spe_targeted <- spe_targeted[, new_order]



## Metrics QC
metrics_qc <- function(spe, spename) {
    qcstats <- perCellQCMetrics(spe, subsets = list(
        Mito = which(seqnames(spe) == "chrM")
    ))

    qc_df <- data.frame(
        log2sum = log2(qcstats$sum),
        log2detected = log2(qcstats$detected),
        subsets_Mito_percent = qcstats$subsets_Mito_percent,
        sample_id = spe$sample_id_short
    )

    qcfilter <- DataFrame(
        low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = spe$sample_id_short),
        low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = spe$sample_id_short),
        high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = spe$sample_id_short)
    )
    qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent


    spe$scran_low_lib_size_low_mito <- factor(qcfilter$low_lib_size & qc_df$subsets_Mito_percent < 0.5, levels = c("TRUE", "FALSE"))


    spe$scran_discard <-
        factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
    spe$scran_low_lib_size <-
        factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
    spe$scran_low_n_features <-
        factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
    spe$scran_high_subsets_Mito_percent <-
        factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

    ## Find edge spots
    spots <- data.frame(
        row = spe$array_row,
        col = spe$array_col,
        sample_id = spe$sample_id_short
    )

    edge_spots_row <- group_by(spots, sample_id, row) %>% summarize(min_col = min(col), max_col = max(col))
    edge_spots_col <- group_by(spots, sample_id, col) %>% summarize(min_row = min(row), max_row = max(row))

    spots <- left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
    spots$edge_spots <- with(spots, row == min_row | row == max_row | col == min_col | col == max_col)

    spots$row_distance <- with(spots, pmin(abs(row - min_row), abs(row - max_row)))
    spots$col_distance <- with(spots, pmin(abs(col - min_col), abs(col - max_col)))
    ## spots$edge_distance <- with(spots, sqrt(row_distance^2 + col_distance^2))
    ## The above is from:
    ## sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
    ## but it was wrong, here's a case the the smallest distance is on the column:
    ## sqrt(0^2 + col_distance^2) = col_distance
    spots$edge_distance <- with(spots, pmin(row_distance, col_distance))


    spe$edge_spots <- factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
    spe$edge_distance <- spots$edge_distance


    spe$scran_low_lib_size_edge <- factor(qcfilter$low_lib_size & spots$edge_distance < 1, levels = c("TRUE", "FALSE"))

    return(spe)
}


spe_wholegenome <- metrics_qc(spe_wholegenome, "wholegenome")
spe_targeted <- metrics_qc(spe_targeted, "targeted")


## Read information about the glare spots manually identified by Madhavi Tippani
## and Sang Ho Kwon.
glare <- read.csv(here("code", "07_spot_qc", "glare_spots.csv"))
colnames(glare)[1] <- "slide"
m <- sapply(
    with(glare, paste0(slide, "_", array)),
    function(partial_sample_id) {
        grep(partial_sample_id, spe_wholegenome$sample_id)[1]
    }
)
glare$sample_id <- spe_wholegenome$sample_id[m]
glare$key <- with(glare, paste0(spotID, "_", sample_id))

## Locate the glare spots
m <- match(glare$key, spe_wholegenome$key)
stopifnot(all(!is.na(m)))
stopifnot(identical(m, match(glare$key, spe_targeted$key)))
# spe_wholegenome <- spe_wholegenome[, -m]
# spe_targeted <- spe_targeted[, -m]
spe_targeted$glare <- spe_wholegenome$glare <- FALSE
spe_targeted$glare[m] <- spe_wholegenome$glare[m] <- TRUE


## Locate low library size spots on the edge for either whole genome or
## targeted sequencing
addmargins(table("wholegenome" = spe_wholegenome$scran_low_lib_size_edge, "targeted" = spe_targeted$scran_low_lib_size_edge))
#            targeted
# wholegenome  TRUE FALSE   Sum
#       TRUE    125    19   144
#       FALSE     8 38115 38123
#       Sum     133 38134 38267

drop_low_library_edge_either <- spe_wholegenome$scran_low_lib_size_edge == "TRUE" | spe_targeted$scran_low_lib_size_edge == "TRUE"
# spe_wholegenome <- spe_wholegenome[, !drop_low_library_edge_either]
# spe_targeted <- spe_targeted[, !drop_low_library_edge_either]

## Clean up some variables names
# spe_targeted$scran_low_lib_size_edge <- spe_wholegenome$scran_low_lib_size_edge <- NULL

spe_wholegenome$drop_low_library_edge_either <- spe_targeted$drop_low_library_edge_either <- drop_low_library_edge_either

table(spe_wholegenome$drop_low_library_edge_either, spe_wholegenome$glare)
#       FALSE  TRUE
# FALSE 38115    20
# TRUE    152     0

table("lq edge" = spe_wholegenome$drop_low_library_edge_either, "glare" = spe_wholegenome$glare, "discard" = spe_wholegenome$scran_discard)
# , , discard = TRUE
#
#        glare
# lq edge FALSE  TRUE
#   FALSE   930     0
#   TRUE    148     0
#
# , , discard = FALSE
#
#        glare
# lq edge FALSE  TRUE
#   FALSE 37185    20
#   TRUE      4     0


spe_wholegenome$quality_groups <- "Pass"
spe_wholegenome$quality_groups[spe_wholegenome$scran_discard == "TRUE"] <- "LQ: retained"
spe_wholegenome$quality_groups[spe_wholegenome$glare] <- "LQ: glare"
spe_wholegenome$quality_groups[spe_wholegenome$drop_low_library_edge_either] <- "LQ: low lib size & edge"
table(spe_wholegenome$quality_groups)
# LQ: glare LQ: low lib size & edge            LQ: retained                    Pass
#        20                     152                     930                   37185

quality_groups_colors <- c("Pass" = "grey90", "LQ: retained" = "orange", "LQ: glare" = "steelblue3", "LQ: low lib size & edge" = "violetred")
p_list <- vis_grid_clus(
    spe = spe_wholegenome,
    clustervar = "quality_groups",
    sort_clust = FALSE,
    colors = quality_groups_colors,
    spatial = FALSE,
    point_size = 2,
    return_plots = TRUE
)

pdf(file.path(dir_plots, "wholegenome_quality_groups.pdf"), useDingbats = FALSE, height = 8 * 4, width = 9 * 3)
print(cowplot::plot_grid(plotlist = p_list, ncol = 3, align = "hv"))
dev.off()

pdf(
    file.path(dir_plots, "wholegenome_quality_groups_library_size.pdf"),
    useDingbats = FALSE,
    width = 18
)
plotColData(spe_wholegenome,
    y = "sum_umi",
    x = "sample_id",
    color_by = "quality_groups") +
    scale_y_log10() +
    scale_colour_manual(values = quality_groups_colors) +
    ylab("Library size (in log10)") +
    theme(axis.title.y = element_text(size = 18))
dev.off()

pdf(file.path(dir_plots, "wholegenome_quality_groups_chrM_proportion.pdf"), useDingbats = FALSE, width = 18)
plotColData(spe_wholegenome,
    y = "expr_chrM_ratio",
    x = "sample_id",
    color_by = "quality_groups") +
    scale_colour_manual(values = quality_groups_colors) +
    ylab("chrM proportion") +
    theme(axis.title.y = element_text(size = 18))
dev.off()

pdf(file.path(dir_plots, "wholegenome_quality_groups_genes_detected.pdf"), useDingbats = FALSE, width = 18)
plotColData(spe_wholegenome,
    y = "sum_gene",
    x = "sample_id",
    color_by = "quality_groups") +
    scale_y_log10() +
    scale_colour_manual(values = quality_groups_colors) +
    ylab("Genes detected (in log10)") +
    theme(axis.title.y = element_text(size = 18))
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 (2023-06-16)
#  os       macOS Ventura 13.5
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2023-08-22
#  rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
#  pandoc   3.1.5 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  abind                    1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
#  AnnotationDbi            1.62.2    2023-07-02 [1] Bioconductor
#  AnnotationHub            3.8.0     2023-04-25 [1] Bioconductor
#  attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.3.0)
#  backports                1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
#  beachmat                 2.16.0    2023-04-25 [1] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.3.0)
#  benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.3.0)
#  benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.3.0)
#  Biobase                * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocFileCache            2.8.0     2023-04-25 [1] Bioconductor
#  BiocGenerics           * 0.46.0    2023-04-25 [1] Bioconductor
#  BiocIO                   1.10.0    2023-04-25 [1] Bioconductor
#  BiocManager              1.30.22   2023-08-08 [1] CRAN (R 4.3.0)
#  BiocNeighbors            1.18.0    2023-04-25 [1] Bioconductor
#  BiocParallel             1.34.2    2023-05-28 [1] Bioconductor
#  BiocSingular             1.16.0    2023-04-25 [1] Bioconductor
#  BiocVersion              3.17.1    2022-12-20 [1] Bioconductor
#  Biostrings               2.68.1    2023-05-16 [1] Bioconductor
#  bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  bluster                  1.10.0    2023-05-08 [1] Bioconductor
#  brio                     1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  broom                    1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
#  bslib                    0.5.1     2023-08-11 [1] CRAN (R 4.3.0)
#  cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr                    3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  car                      3.1-2     2023-03-30 [1] CRAN (R 4.3.0)
#  carData                  3.0-5     2022-01-06 [1] CRAN (R 4.3.0)
#  cli                      3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  cluster                  2.1.4     2022-08-22 [1] CRAN (R 4.3.1)
#  codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
#  colorout                 1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
#  config                   0.3.1     2020-12-17 [1] CRAN (R 4.3.0)
#  cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
#  crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  curl                     5.0.2     2023-08-14 [1] CRAN (R 4.3.1)
#  data.table               1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
#  DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
#  dbplyr                   2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
#  DelayedArray             0.26.7    2023-07-30 [1] Bioconductor
#  DelayedMatrixStats       1.22.5    2023-08-10 [1] Bioconductor
#  devtools               * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
#  digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
#  doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.3.0)
#  dotCall64                1.0-2     2022-10-03 [1] CRAN (R 4.3.0)
#  dplyr                  * 1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
#  dqrng                    0.3.0     2021-05-01 [1] CRAN (R 4.3.0)
#  DropletUtils             1.20.0    2023-05-08 [1] Bioconductor
#  DT                       0.28      2023-05-18 [1] CRAN (R 4.3.0)
#  edgeR                    3.42.4    2023-05-31 [1] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  ExperimentHub            2.8.1     2023-07-16 [1] Bioconductor
#  fansi                    1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
#  farver                   2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
#  fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  fields                   14.1      2022-08-12 [1] CRAN (R 4.3.0)
#  filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
#  foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
#  fs                       1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
#  generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb           * 1.36.1    2023-07-02 [1] Bioconductor
#  GenomeInfoDbData         1.2.10    2023-05-06 [1] Bioconductor
#  GenomicAlignments        1.36.0    2023-04-25 [1] Bioconductor
#  GenomicRanges          * 1.52.0    2023-04-25 [1] Bioconductor
#  ggbeeswarm               0.7.2     2023-04-29 [1] CRAN (R 4.3.0)
#  ggplot2                * 3.4.3     2023-08-14 [1] CRAN (R 4.3.1)
#  ggpubr                 * 0.6.0     2023-02-10 [1] CRAN (R 4.3.0)
#  ggrepel                  0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
#  ggsignif                 0.6.4     2022-10-13 [1] CRAN (R 4.3.0)
#  glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  golem                    0.4.1     2023-06-05 [1] CRAN (R 4.3.0)
#  gridExtra                2.3       2017-09-09 [1] CRAN (R 4.3.0)
#  gtable                   0.3.3     2023-03-21 [1] CRAN (R 4.3.0)
#  HDF5Array                1.28.1    2023-05-01 [1] Bioconductor
#  here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                      1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools                0.5.6     2023-08-10 [1] CRAN (R 4.3.0)
#  htmlwidgets              1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv                   1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
#  httr                     1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
#  igraph                   1.5.1     2023-08-10 [1] CRAN (R 4.3.0)
#  interactiveDisplayBase   1.38.0    2023-04-25 [1] Bioconductor
#  IRanges                * 2.34.1    2023-07-02 [1] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
#  iterators                1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
#  jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.3.0)
#  jsonlite                 1.8.7     2023-06-29 [1] CRAN (R 4.3.0)
#  KEGGREST                 1.40.0    2023-04-25 [1] Bioconductor
#  labeling                 0.4.2     2020-10-20 [1] CRAN (R 4.3.0)
#  later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
#  lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.3.1)
#  lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
#  lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
#  limma                    3.56.2    2023-06-04 [1] Bioconductor
#  locfit                   1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
#  lubridate                1.9.2     2023-02-10 [1] CRAN (R 4.3.0)
#  magick                   2.7.5     2023-08-07 [1] CRAN (R 4.3.0)
#  magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  maps                     3.4.1     2022-10-30 [1] CRAN (R 4.3.0)
#  Matrix                   1.6-1     2023-08-14 [1] CRAN (R 4.3.0)
#  MatrixGenerics         * 1.12.3    2023-07-30 [1] Bioconductor
#  matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
#  memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  metapod                  1.8.0     2023-04-25 [1] Bioconductor
#  mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
#  paletteer                1.5.0     2022-10-19 [1] CRAN (R 4.3.0)
#  pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild                 1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
#  pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                  1.3.2.1   2023-07-08 [1] CRAN (R 4.3.0)
#  plotly                   4.10.2    2023-06-03 [1] CRAN (R 4.3.0)
#  png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
#  prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
#  processx                 3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
#  profvis                  0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  promises                 1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                   1.0.1     2023-05-06 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                       1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
#  purrr                    1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
#  R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
#  R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
#  R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
#  R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.3.0)
#  remotes                  2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
#  restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
#  rhdf5                    2.44.0    2023-04-25 [1] Bioconductor
#  rhdf5filters             1.12.1    2023-04-30 [1] Bioconductor
#  Rhdf5lib                 1.22.0    2023-04-25 [1] Bioconductor
#  rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
#  rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
#  Rsamtools                2.16.0    2023-04-25 [1] Bioconductor
#  RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
#  rstatix                  0.7.2     2023-02-01 [1] CRAN (R 4.3.0)
#  rsthemes                 0.4.0     2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
#  rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
#  rtracklayer              1.60.0    2023-04-25 [1] Bioconductor
#  S4Arrays                 1.0.5     2023-07-24 [1] Bioconductor
#  S4Vectors              * 0.38.1    2023-05-02 [1] Bioconductor
#  sass                     0.4.7     2023-07-15 [1] CRAN (R 4.3.0)
#  ScaledMatrix             1.8.1     2023-05-03 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
#  scater                 * 1.28.0    2023-04-25 [1] Bioconductor
#  scran                  * 1.28.2    2023-07-23 [1] Bioconductor
#  scuttle                * 1.10.2    2023-08-03 [1] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                    1.7.5     2023-08-12 [1] CRAN (R 4.3.0)
#  shinyWidgets             0.7.6     2023-01-08 [1] CRAN (R 4.3.0)
#  SingleCellExperiment   * 1.22.0    2023-04-25 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [1] CRAN (R 4.3.0)
#  sparseMatrixStats        1.12.2    2023-07-02 [1] Bioconductor
#  SpatialExperiment      * 1.10.0    2023-04-25 [1] Bioconductor
#  spatialLIBD            * 1.13.4    2023-05-24 [1] Github (LieberInstitute/spatialLIBD@edc8b72)
#  statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
#  stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
#  stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
#  SummarizedExperiment   * 1.30.2    2023-06-11 [1] Bioconductor
#  suncalc                  0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  testthat               * 3.1.10    2023-07-06 [1] CRAN (R 4.3.0)
#  tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
#  tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  timechange               0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
#  urlchecker               1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis                * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
#  utf8                     1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs                    0.6.3     2023-06-14 [1] CRAN (R 4.3.0)
#  vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.3.0)
#  viridis                  0.6.4     2023-07-22 [1] CRAN (R 4.3.0)
#  viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
#  withr                    2.5.0     2022-03-03 [1] CRAN (R 4.3.0)
#  XML                      3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
#  xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                  0.40.0    2023-04-25 [1] Bioconductor
#  yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
#  zlibbioc                 1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
