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

p_list <- vis_grid_clus(
    spe = spe_wholegenome,
    clustervar = "quality_groups",
    sort_clust = FALSE,
    colors = c("Pass" = "grey90", "LQ: retained" = "orange", "LQ: glare" = "steelblue3", "LQ: low lib size & edge" = "violetred"),
    spatial = FALSE,
    point_size = 2,
    return_plots = TRUE
)

pdf(file.path(dir_plots, "wholegenome_quality_groups.pdf"), useDingbats = FALSE, height = 8 * 4, width = 9 * 3)
print(cowplot::plot_grid(plotlist = p_list, ncol = 3, align = "hv"))
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
