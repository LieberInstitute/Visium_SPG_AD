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
load(here::here("processed-data", "spe", "spe.Rdata"), verbose = TRUE)
load(here::here("processed-data", "spe", "spe_targeted.Rdata"), verbose = TRUE)

## Create output directories
dir_plots <- here::here("plots", "07_spot_qc")
dir_rdata <- here::here("processed-data", "07_spot_qc")
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

slide_order <- c("V10A27106", "V10T31036", "V10A27004")
sample_order <- unlist(sapply(slide_order, function(i) {
    sort(unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
sample_order

shorten_names <- function(x) {
    gsub(slide_order[3], "S3", gsub(slide_order[2], "S2", gsub(slide_order[1], "S1", x)))
}

sample_id_factor <- factor(shorten_names(spe$sample_id), levels = shorten_names(sample_order))

palette_colors <- paste0(
    rep(c("darkorchid", "sienna", "steelblue"), c(4, 4, 2)),
    c(1:4, 1:4, 3:4)
)
names(palette_colors) <- levels(sample_id_factor)

## Metrics QC
metrics_qc <- function(spe, spename) {
    qcstats <- perCellQCMetrics(spe, subsets = list(
        Mito = which(seqnames(spe) == "chrM")
    ))

    qc_df <- data.frame(
        log2sum = log2(qcstats$sum),
        log2detected = log2(qcstats$detected),
        subsets_Mito_percent = qcstats$subsets_Mito_percent,
        sample_id = sample_id_factor
    )

    for (i in c("log2sum", "log2detected", "subsets_Mito_percent")) {
        pdf(file.path(dir_plots, paste0("scran_", spename, "_distribution_", i, ".pdf")), width = 14)
        p <- ggpubr::ggviolin(
            qc_df,
            "sample_id",
            i,
            fill = "sample_id",
            palette = palette_colors,
            add = "boxplot",
            add.params = list(fill = "white")
        )
        print(p)
        dev.off()
    }

    qcfilter <- DataFrame(
        low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = sample_id_factor),
        low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = sample_id_factor),
        high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = sample_id_factor)
    )
    qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent


    pdf(file.path(
        dir_plots,
        paste0("scran_", spename, "_low_lib_size_vs_mito_scatter.pdf")
    ), width = 21)
    ggplot(
        qc_df,
        aes(
            x = log2sum,
            y = subsets_Mito_percent,
            color = qcfilter$low_lib_size,
            shape = isOutlier(qcstats$subsets_Mito_percent, type = "lower", batch = sample_id_factor)
        )
    ) +
        geom_point() +
        facet_grid(~sample_id) +
        guides(color = guide_legend("Low lib?")) +
        guides(shape = guide_legend("Low mito?"))
    dev.off()

    spe$scran_low_lib_size_low_mito <- factor(qcfilter$low_lib_size & qc_df$subsets_Mito_percent < 0.5, levels = c("TRUE", "FALSE"))
    vis_grid_clus(
        spe = spe,
        clustervar = "scran_low_lib_size_low_mito",
        pdf = file.path(dir_plots, paste0("scran_", spename, "_low_lib_size_vs_mito.pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        spatial = FALSE,
        point_size = 2,
        sample_order = sample_order
    )


    spe$scran_discard <-
        factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
    spe$scran_low_lib_size <-
        factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
    spe$scran_low_n_features <-
        factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
    spe$scran_high_subsets_Mito_percent <-
        factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

    for (i in colnames(qcfilter)) {
        vis_grid_clus(
            spe = spe,
            clustervar = paste0("scran_", i),
            pdf = file.path(dir_plots, paste0("scran_", spename, "_", i, ".pdf")),
            sort_clust = FALSE,
            colors = c("FALSE" = "grey90", "TRUE" = "orange"),
            spatial = FALSE,
            point_size = 2,
            sample_order = sample_order
        )
    }


    ## Find edge spots
    spots <- data.frame(
        row = spatialData(spe)$array_row,
        col = spatialData(spe)$array_col,
        sample_id = sample_id_factor
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
    vis_grid_clus(
        spe = spe,
        clustervar = "edge_spots",
        pdf = file.path(dir_plots, paste0("egde_spots_", spename, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        spatial = FALSE,
        point_size = 2,
        sample_order = sample_order
    )

    vis_grid_gene(
        spe = spe,
        geneid = "edge_distance",
        pdf = file.path(dir_plots, paste0("egde_distance_", spename, ".pdf")),
        spatial = FALSE,
        point_size = 2,
        sample_order = sample_order
    )


    pdf(file.path(
        dir_plots,
        paste0("scran_", spename, "_low_lib_size_vs_edge_distance_scatter.pdf")
    ), width = 21)
    ggplot(
        qc_df,
        aes(
            x = log2sum,
            y = spots$edge_distance,
            color = qcfilter$low_lib_size
        )
    ) +
        geom_point() +
        facet_grid(~sample_id) +
        guides(color = guide_legend("Low lib?"))
    dev.off()

    spe$scran_low_lib_size_edge <- factor(qcfilter$low_lib_size & spots$edge_distance < 1, levels = c("TRUE", "FALSE"))
    vis_grid_clus(
        spe = spe,
        clustervar = "scran_low_lib_size_edge",
        pdf = file.path(dir_plots, paste0("scran_", spename, "_low_lib_size_vs_edge_distance.pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        spatial = FALSE,
        point_size = 2,
        sample_order = sample_order
    )

    return(spe)
}


spe <- metrics_qc(spe, "wholegenome")
spe_targeted <- metrics_qc(spe_targeted, "targeted")


## Segmentation spot QC
seg_df <- data.frame(
    Percent_Abeta = spe$PAbeta,
    Percent_DAPI = spe$PDAPI,
    Percent_pTau = spe$PpTau,
    Number_Abeta = spe$NAbeta,
    Number_DAPI = spe$NDAPI,
    Number_pTau = spe$NpTau,
    sample_id = sample_id_factor,
    check.names = FALSE
)

for (i in paste0(rep(c("Percent_", "Number_"), each = 3), rep(c("Abeta", "DAPI", "pTau"), 2))) {
    pdf(file.path(dir_plots, paste0("segmentation_distribution_", i, ".pdf")), width = 14)
    p <- ggpubr::ggviolin(
        seg_df,
        "sample_id",
        i,
        fill = "sample_id",
        palette = palette_colors,
        add = "boxplot",
        add.params = list(fill = "white")
    )
    print(p)
    dev.off()
}

segqcfilter <- DataFrame(
    Number_DAPI = isOutlier(spe$NDAPI, type = "higher", batch = sample_id_factor, nmads = 3),
    Percent_pTau = isOutlier(spe$PpTau, type = "higher", batch = sample_id_factor, nmads = 5)
)
segqcfilter$discard_segmentation <- segqcfilter$Number_DAPI | segqcfilter$Percent_pTau

attr(segqcfilter$Number_DAPI, "thresholds")
#        S1_A1_Br3874 S1_B1_Br3854 S1_C1_Br3873 S1_D1_Br3880 S2_A1_Br3874
# lower          -Inf         -Inf         -Inf         -Inf         -Inf
# higher      11.8956      12.8956      12.8956      13.8956      12.8956
#        S2_B1_Br3854 S2_C1_Br3873 S2_D1_Br3880 S3_A1_Br3874 S3_D1_Br3880
# lower          -Inf         -Inf         -Inf         -Inf         -Inf
# higher      12.8956      12.8956      13.8956       7.4478      14.8956
attr(segqcfilter$Percent_pTau, "thresholds")
#        S1_A1_Br3874 S1_B1_Br3854 S1_C1_Br3873 S1_D1_Br3880 S2_A1_Br3874
# lower          -Inf         -Inf         -Inf         -Inf         -Inf
# higher            0   0.01437132            0   0.02798106    0.0012306
#        S2_B1_Br3854 S2_C1_Br3873 S2_D1_Br3880 S3_A1_Br3874 S3_D1_Br3880
# lower          -Inf         -Inf         -Inf         -Inf         -Inf
# higher   0.03801507            0   0.04665947            0   0.03044226

spe$scran_discard_segmentation <-
    factor(segqcfilter$discard_segmentation, levels = c("TRUE", "FALSE"))
spe$scran_Number_DAPI <-
    factor(segqcfilter$Number_DAPI, levels = c("TRUE", "FALSE"))
spe$scran_Percent_pTau <-
    factor(segqcfilter$Percent_pTau, levels = c("TRUE", "FALSE"))

for (i in colnames(segqcfilter)) {
    vis_grid_clus(
        spe = spe,
        clustervar = paste0("scran_", i),
        pdf = file.path(dir_plots, paste0("scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        spatial = FALSE,
        point_size = 2
    )
}

## Don't drop any spots due to isOutlier() results on the image segmentation
## Drop this information since we won't use for downstream analyses.
spe$scran_discard_segmentation <- NULL
spe$scran_Number_DAPI <- NULL
spe$scran_Percent_pTau <- NULL

## Read information about the glare spots manually identified by Madhavi Tippani
## and Sang Ho Kwon.
glare <- read.csv(here("code", "07_spot_qc", "glare_spots.csv"))
colnames(glare)[1] <- "slide"
m <- sapply(
    with(glare, paste0(slide, "_", array)),
    function(partial_sample_id) {
        grep(partial_sample_id, spe$sample_id)[1]
    }
)
glare$sample_id <- spe$sample_id[m]
glare$key <- with(glare, paste0(spotID, "_", sample_id))

## Locate and drop the glare spots
m <- match(glare$key, spe$key)
stopifnot(all(!is.na(m)))
stopifnot(identical(m, match(glare$key, spe_targeted$key)))
spe <- spe[, -m]
spe_targeted <- spe_targeted[, -m]


## Drop low library size spots on the edge for either whole genome or
## targeted sequencing
addmargins(table("wholegenome" = spe$scran_low_lib_size_edge, "targeted" = spe_targeted$scran_low_lib_size_edge))
#            targeted
# wholegenome  TRUE FALSE   Sum
#       TRUE    125    19   144
#       FALSE     8 38115 38123
#       Sum     133 38134 38267

drop_low_library_edge_either <- spe$scran_low_lib_size_edge == "TRUE" | spe_targeted$scran_low_lib_size_edge == "TRUE"
spe <- spe[, !drop_low_library_edge_either]
spe_targeted <- spe_targeted[, !drop_low_library_edge_either]

## Clean up some variables names
spe_targeted$scran_low_lib_size_edge <- spe$scran_low_lib_size_edge <- NULL

## Save for later
save(spe, file = here::here("processed-data", "spe", "spe_postqc.Rdata"))
save(spe_targeted,
    file = here::here("processed-data", "spe", "spe_targeted_postqc.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.2 Patched (2021-11-04 r81138)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-01-12
#  pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date (UTC) lib source
#  abind                    1.4-5    2016-07-21 [2] CRAN (R 4.1.0)
#  AnnotationDbi            1.56.2   2021-11-09 [2] Bioconductor
#  AnnotationHub            3.2.0    2021-10-26 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.2)
#  backports                1.4.1    2021-12-13 [2] CRAN (R 4.1.2)
#  beachmat                 2.10.0   2021-10-26 [2] Bioconductor
#  beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.2)
#  benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.2)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.2)
#  Biobase                * 2.54.0   2021-10-26 [2] Bioconductor
#  BiocFileCache            2.2.0    2021-10-26 [2] Bioconductor
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
#  bluster                  1.4.0    2021-10-26 [1] Bioconductor
#  broom                    0.7.11   2022-01-03 [2] CRAN (R 4.1.2)
#  bslib                    0.3.1    2021-10-06 [2] CRAN (R 4.1.2)
#  cachem                   1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
#  callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
#  car                      3.0-12   2021-11-06 [2] CRAN (R 4.1.2)
#  carData                  3.0-5    2022-01-06 [2] CRAN (R 4.1.2)
#  cli                      3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
#  cluster                  2.1.2    2021-04-17 [3] CRAN (R 4.1.2)
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
#  dplyr                  * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.2)
#  DropletUtils             1.14.1   2021-11-08 [1] Bioconductor
#  DT                       0.20     2021-11-15 [2] CRAN (R 4.1.2)
#  edgeR                    3.36.0   2021-10-26 [2] Bioconductor
#  ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  ExperimentHub            2.2.0    2021-10-26 [2] Bioconductor
#  fansi                    1.0.0    2022-01-10 [2] CRAN (R 4.1.2)
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   13.3     2021-10-30 [2] CRAN (R 4.1.2)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                       1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
#  generics                 0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
#  GenomeInfoDb           * 1.30.0   2021-10-26 [2] Bioconductor
#  GenomeInfoDbData         1.2.7    2021-11-01 [2] Bioconductor
#  GenomicAlignments        1.30.0   2021-10-26 [2] Bioconductor
#  GenomicRanges          * 1.46.1   2021-11-18 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.2)
#  ggplot2                * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  ggpubr                 * 0.4.0    2020-06-27 [2] CRAN (R 4.1.2)
#  ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  ggsignif                 0.6.3    2021-09-09 [2] CRAN (R 4.1.2)
#  glue                     1.6.0    2021-12-17 [2] CRAN (R 4.1.2)
#  golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.2)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.22.1   2021-11-14 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.2)
#  htmltools                0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
#  htmlwidgets              1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
#  httpuv                   1.6.5    2022-01-05 [2] CRAN (R 4.1.2)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  igraph                   1.2.11   2022-01-04 [2] CRAN (R 4.1.2)
#  interactiveDisplayBase   1.32.0   2021-10-26 [2] Bioconductor
#  IRanges                * 2.28.0   2021-10-26 [2] Bioconductor
#  irlba                    2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  KEGGREST                 1.34.0   2021-10-26 [2] Bioconductor
#  knitr                    1.37     2021-12-16 [2] CRAN (R 4.1.2)
#  later                    1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
#  lattice                  0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
#  limma                    3.50.0   2021-10-26 [2] Bioconductor
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
#  magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  maps                     3.4.0    2021-09-25 [2] CRAN (R 4.1.2)
#  Matrix                   1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
#  MatrixGenerics         * 1.6.0    2021-10-26 [2] Bioconductor
#  matrixStats            * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
#  memoise                  2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
#  metapod                  1.2.0    2021-10-26 [1] Bioconductor
#  mime                     0.12     2021-09-28 [2] CRAN (R 4.1.2)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
#  pkgbuild                 1.3.1    2021-12-20 [2] CRAN (R 4.1.2)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.2.4    2021-11-30 [2] CRAN (R 4.1.2)
#  plotly                   4.10.0   2021-10-09 [2] CRAN (R 4.1.2)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome               1.3.1    2021-07-16 [1] CRAN (R 4.1.2)
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
#  Rcpp                     1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                    1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
#  remotes                  2.4.2    2021-11-30 [2] CRAN (R 4.1.2)
#  restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                    2.38.0   2021-10-26 [2] Bioconductor
#  rhdf5filters             1.6.0    2021-10-26 [2] Bioconductor
#  Rhdf5lib                 1.16.0   2021-10-26 [2] Bioconductor
#  rjson                    0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
#  rlang                    0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
#  rmote                    0.3.4    2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.1.2    2021-09-08 [2] CRAN (R 4.1.2)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools                2.10.0   2021-10-26 [2] Bioconductor
#  RSQLite                  2.2.9    2021-12-06 [2] CRAN (R 4.1.2)
#  rstatix                  0.7.0    2021-02-13 [2] CRAN (R 4.1.2)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.2)
#  rtracklayer              1.54.0   2021-10-26 [2] Bioconductor
#  S4Vectors              * 0.32.3   2021-11-21 [2] Bioconductor
#  sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
#  ScaledMatrix             1.2.0    2021-10-26 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                 * 1.22.0   2021-10-26 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.2)
#  scran                  * 1.22.1   2021-11-14 [1] Bioconductor
#  scuttle                * 1.4.0    2021-10-26 [1] Bioconductor
#  servr                    0.24     2021-11-16 [1] CRAN (R 4.1.2)
#  sessioninfo            * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
#  shiny                    1.7.1    2021-10-02 [2] CRAN (R 4.1.2)
#  shinyWidgets             0.6.2    2021-09-17 [1] CRAN (R 4.1.2)
#  SingleCellExperiment   * 1.16.0   2021-10-26 [2] Bioconductor
#  spam                     2.8-0    2022-01-06 [2] CRAN (R 4.1.2)
#  sparseMatrixStats        1.6.0    2021-10-26 [2] Bioconductor
#  SpatialExperiment      * 1.4.0    2021-10-26 [1] Bioconductor
#  spatialLIBD            * 1.6.4    2021-12-08 [1] Github (LieberInstitute/spatialLIBD@ea2c037)
#  statmod                  1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
#  stringi                  1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.24.0   2021-10-26 [2] Bioconductor
#  testthat                 3.1.1    2021-12-03 [2] CRAN (R 4.1.2)
#  tibble                   3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
#  tidyr                    1.1.4    2021-09-27 [2] CRAN (R 4.1.2)
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
#  yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc                 1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
