## import required libraries
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
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)
load(here::here("processed-data", "07_spot_qc", "spe_targeted_postqc.Rdata"), verbose = TRUE)



## output directories
dir_plots <- here::here("plots", "07_spot_qc", "outliers")
dir.create(dir_plots, showWarnings = FALSE)


## Create scatter plots to view outlier information
create_plots <- function(spe_object, pathology, n = 0, p = 0.01) {
    # pathology can be 'Abeta' or 'pTau'
    # n = threshold for number of pathology 'blobs' in spot
    # p = percentage of pathology pixels in spot
    # add optional params (diagnosis)

    path_df <- data.frame(
        spot_id = rownames(colData(spe_object)),
        diagnosis = colData(spe_object)$diagnosis,
        sample_id = colData(spe_object)$sample_id,
        NAbeta = colData(spe_object)$NAbeta,
        NpTau = colData(spe_object)$NpTau,
        PAbeta = colData(spe_object)$PAbeta,
        PpTau = colData(spe_object)$PpTau
    )

    ## row 1: Percent > 0.01, column 2: Percent absent
    if (pathology == "Abeta") {
        path_df <- path_df |> mutate(outliers = case_when(
            NAbeta > n & PAbeta > p ~ "both",
            NAbeta > n & PAbeta <= p ~ "n",
            NAbeta <= n & PAbeta > p ~ "p",
            NAbeta <= n & PAbeta <= p ~ "none"
        ))

        plot <- ggpubr::ggscatter(path_df,
            x = "NAbeta", y = "PAbeta",
            color = "outliers", size = 0.5
        )
        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )
    }

    if (pathology == "pTau") {
        path_df <- path_df |> mutate(outliers = case_when(
            NpTau > n & PpTau > p ~ "both",
            NpTau > n & PpTau <= p ~ "n",
            NpTau <= n & PpTau > p ~ "p",
            NpTau <= n & PpTau <= p ~ "none"
        ))

        plot <- ggpubr::ggscatter(path_df,
            x = "NpTau", y = "PpTau",
            color = "outliers", size = 0.5
        )
        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )
    }

    return(plot)
}



# create plots for whole genome
for (path in c("Abeta", "pTau")) {
    pdf(file.path(dir_plots, paste0("spe_whole", "_", path, "_", "scatter", ".pdf")), width = 14)
    print(create_plots(spe, path))
    dev.off()
}

# create plots for targeted genome
for (path in c("Abeta", "pTau")) {
    pdf(file.path(dir_plots, paste0("spe_targeted", "_", path, "_", "scatter", ".pdf")), width = 14)
    print(create_plots(spe_targeted, path))
    dev.off()
}




## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
