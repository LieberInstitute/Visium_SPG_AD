# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "03_scatter_plots_case_vs_controls",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "10G")

library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("getopt")

## Specify parameters
spec <- matrix(c(
    "speopt$spetype", "s", 2, "character", "SPE opt$spetype: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}



## output plots directory

# dir.create(here::here("plots","07_spot_qc", "outliers" type), showWarnings = FALSE)
dir_plots <- here::here("plots", "07_spot_qc", "outliers")



## Create scatter plots to view outlier information
create_plots <- function(spe_object, pathology) {
    # pathology can be 'Abeta' or 'pTau'
    # n = threshold for number of pathology 'blobs' in spot
    # p = percentage of pathology pixels in spot
    # add optional params (diagnosis)
    colors_hex <- c("#00AFBB", "#E7B800", "#FC4E07", "#CC79A7")

    path_df <- data.frame(
        spot_id = rownames(colData(spe_object)),
        diagnosis = colData(spe_object)$diagnosis,
        sample_id = colData(spe_object)$sample_id,
        NAbeta = as.integer(colData(spe_object)$NAbeta),
        NpTau = as.integer(colData(spe_object)$NpTau),
        PAbeta = colData(spe_object)$PAbeta,
        PpTau = colData(spe_object)$PpTau
    )

    path_df$diagnosis <- factor(path_df$diagnosis, levels = c("Control", "AD"))
    ## row 1: Percent > 0.01, column 2: Percent absent
    if (pathology == "Abeta") {
        n = 1
        p = 0.108
        path_df <- path_df |> mutate(outliers = case_when(
            NAbeta > n & PAbeta > p ~ "both",
            NAbeta > n & PAbeta <= p ~ "n",
            NAbeta <= n & PAbeta > p ~ "p",
            NAbeta <= n & PAbeta <= p ~ "none"
        ))

        plot <- ggpubr::ggscatter(path_df,
            x = "NAbeta", y = "PAbeta",
            color = "outliers", size = 0.5, alpha = 0.3
            # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#CC79A7")
        )
        plot <- facet(plot + theme_bw(),
            facet.by = "diagnosis",
            short.panel.labs = TRUE
        )
        plot <- plot + scale_color_manual(
            name = "type of Abeta",
            labels = c("n and %", "n", "p" ,"none"),
            values = colors_hex
        )
    }

    if (pathology == "pTau") {
        n= 8
        p = 0.0143
        path_df <- path_df |> mutate(outliers = case_when(
            NpTau > n & PpTau > p ~ "both",
            NpTau > n & PpTau <= p ~ "n",
            NpTau <= n & PpTau > p ~ "p",
            NpTau <= n & PpTau <= p ~ "none"
        ))

        plot <- ggpubr::ggscatter(path_df,
            x = "NpTau", y = "PpTau",
            color = "outliers", size = 0.5, alpha = 0.3
        )

        plot <- facet(plot + theme_bw(),
            facet.by = "diagnosis",
            short.panel.labs = TRUE
        )

        plot <- plot + scale_color_manual(
            name = "type of pTau",
            labels = c("n and %", "n", "p" ,"none"),
            values = colors_hex
        )
    }

    return(plot)
}



for (type in genome_type) {
    spe <-
        readRDS(
            here::here(
                "processed-data",
                "08_harmony_BayesSpace",
                type,
                paste0("spe_harmony_", type, ".rds")
            )
        )

    spe <- cluster_import(
        spe,
        cluster_dir = here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            type,
            "clusters_BayesSpace"
        ),
        prefix = ""
    )
}


create_plots(spe, "Abeta")
create_plots(spe, "pTau")
