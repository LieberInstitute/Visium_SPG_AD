library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
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

    for(i in c("log2sum", "log2detected", "subsets_Mito_percent")) {
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

for(i in paste0(rep(c("Percent_", "Number_"), each = 3), rep(c("Abeta", "DAPI", "pTau"), 2))) {
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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
