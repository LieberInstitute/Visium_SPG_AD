# sgejobs::job_single(
#     "qc_metrics_and_segmentation",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "80G",
#     command = "Rscript qc_metrics_and_segmentation.R",
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
        low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = spe$sample_id_short),
        low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = spe$sample_id_short),
        high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = spe$sample_id_short)
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
            shape = isOutlier(qcstats$subsets_Mito_percent, type = "lower", batch = spe$sample_id_short)
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
        point_size = 2
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
            point_size = 2
        )
    }


    ## Find edge spots
    spots <- data.frame(
        row = spatialData(spe)$array_row,
        col = spatialData(spe)$array_col,
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
    vis_grid_clus(
        spe = spe,
        clustervar = "edge_spots",
        pdf = file.path(dir_plots, paste0("egde_spots_", spename, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        spatial = FALSE,
        point_size = 2
    )

    vis_grid_gene(
        spe = spe,
        geneid = "edge_distance",
        pdf = file.path(dir_plots, paste0("egde_distance_", spename, ".pdf")),
        spatial = FALSE,
        point_size = 2
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
        point_size = 2
    )

    return(spe)
}


spe_wholegenome <- metrics_qc(spe_wholegenome, "wholegenome")
spe_targeted <- metrics_qc(spe_targeted, "targeted")


## Segmentation spot QC
seg_df <- data.frame(
    Percent_Abeta = spe_wholegenome$PAbeta,
    Percent_DAPI = spe_wholegenome$PDAPI,
    Percent_pTau = spe_wholegenome$PpTau,
    Number_Abeta = spe_wholegenome$NAbeta,
    Number_DAPI = spe_wholegenome$NDAPI,
    Number_pTau = spe_wholegenome$NpTau,
    sample_id = spe_wholegenome$sample_id_short,
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
    Number_DAPI = isOutlier(spe_wholegenome$NDAPI, type = "higher", batch = spe_wholegenome$sample_id_short, nmads = 3),
    Percent_pTau = isOutlier(spe_wholegenome$PpTau, type = "higher", batch = spe_wholegenome$sample_id_short, nmads = 5)
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

spe_wholegenome$scran_discard_segmentation <-
    factor(segqcfilter$discard_segmentation, levels = c("TRUE", "FALSE"))
spe_wholegenome$scran_Number_DAPI <-
    factor(segqcfilter$Number_DAPI, levels = c("TRUE", "FALSE"))
spe_wholegenome$scran_Percent_pTau <-
    factor(segqcfilter$Percent_pTau, levels = c("TRUE", "FALSE"))

for (i in colnames(segqcfilter)) {
    vis_grid_clus(
        spe = spe_wholegenome,
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
spe_wholegenome$scran_discard_segmentation <- NULL
spe_wholegenome$scran_Number_DAPI <- NULL
spe_wholegenome$scran_Percent_pTau <- NULL

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

## Locate and drop the glare spots
m <- match(glare$key, spe_wholegenome$key)
stopifnot(all(!is.na(m)))
stopifnot(identical(m, match(glare$key, spe_targeted$key)))
spe_wholegenome <- spe_wholegenome[, -m]
spe_targeted <- spe_targeted[, -m]


## Drop low library size spots on the edge for either whole genome or
## targeted sequencing
addmargins(table("wholegenome" = spe_wholegenome$scran_low_lib_size_edge, "targeted" = spe_targeted$scran_low_lib_size_edge))
#            targeted
# wholegenome  TRUE FALSE   Sum
#       TRUE    125    19   144
#       FALSE     8 38115 38123
#       Sum     133 38134 38267

drop_low_library_edge_either <- spe_wholegenome$scran_low_lib_size_edge == "TRUE" | spe_targeted$scran_low_lib_size_edge == "TRUE"
spe_wholegenome <- spe_wholegenome[, !drop_low_library_edge_either]
spe_targeted <- spe_targeted[, !drop_low_library_edge_either]

## Clean up some variables names
spe_targeted$scran_low_lib_size_edge <- spe_wholegenome$scran_low_lib_size_edge <- NULL

## Save for later
saveRDS(spe_wholegenome, file.path(dir_rdata, "spe_wholegenome_postqc.rds"))
saveRDS(spe_targeted, file = file.path(dir_rdata, "spe_targeted_postqc.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
