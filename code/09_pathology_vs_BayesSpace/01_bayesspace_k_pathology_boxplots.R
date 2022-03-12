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
spe_wholegenome <- readRDS(
    here::here(
        "processed-data", "07_spot_qc", "spe_wholegenome_postqc.rds"
    )
)
spe_targeted <- readRDS(
    here::here(
        "processed-data", "07_spot_qc", "spe_targeted_postqc.rds"
    )
)

## output directories
dir_plots <- here::here("plots", "09_pathology_vs_BayesSpace", "pathology_vs_Bayesspace_cluster_boxplots")
dir.create(dir_plots, showWarnings = FALSE)

# import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") # , suffix

spe_wholegenome <- cluster_import(
    spe_wholegenome,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = ""
)

# import cluster info for targeted genome
dir_rdata_targeted <- here::here("processed-data", "08_harmony_BayesSpace", "targeted")

spe_targeted <- cluster_import(
    spe_targeted,
    cluster_dir = file.path(dir_rdata_targeted, "clusters_BayesSpace"),
    prefix = ""
)


# create df with relevant variables for whole
cluster_whole_df <- data.frame(
    spot_id = colnames(spe_wholegenome),
    diagnosis = spe_wholegenome$diagnosis,
    sample_id = factor(spe_wholegenome$sample_id, levels = unique(spe_wholegenome$sample_id)),
    NAbeta = spe_wholegenome$NAbeta,
    NpTau = spe_wholegenome$NpTau,
    PAbeta = spe_wholegenome$PAbeta,
    PpTau = spe_wholegenome$PpTau,
    colData(spe_wholegenome)[grep("BayesSpace_harmony_k", colnames(colData(spe_wholegenome)))] ## could use select(matches =) but
    ## have to convert entire colData to df
)

# create df with relevant variables for targeted
cluster_targeted_df <- data.frame(
    spot_id = colnames(spe_targeted),
    diagnosis = spe_targeted$diagnosis,
    sample_id = factor(spe_targeted$sample_id, levels = unique(spe_targeted$sample_id)),
    NAbeta = spe_targeted$NAbeta,
    NpTau = spe_targeted$NpTau,
    PAbeta = spe_targeted$PAbeta,
    PpTau = spe_targeted$PpTau,
    colData(spe_targeted)[grep("BayesSpace_harmony_k", colnames(colData(spe_targeted)))]
)


# convert cluster info to factor in cluster_whole_df
cols_whole <- colnames(cluster_whole_df |> select(matches("harmony")))
cluster_whole_df[cols_whole] <- lapply(cluster_whole_df[cols_whole], factor)

# convert cluster info to factor in cluster_targeted_df
cols_targeted <- colnames(cluster_targeted_df |> select(matches("harmony")))
cluster_targeted_df[cols_targeted] <- lapply(cluster_targeted_df[cols_targeted], factor)



pathology_measures <- c("PpTau", "PAbeta")

# create plots for whole genome
for (measure in pathology_measures) {
    pdf(file.path(dir_plots, paste0("spe_wholegenome", "_", measure, ".pdf")), width = 18)

    for (i in cols_whole) {
        plot <- ggpubr::ggviolin(
            cluster_whole_df,
            i,
            measure,
            add = "boxplot",
            add.params = list(fill = "white")
        )
        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE,
            scales = "free_y"
        )


        print(plot)
    }
    dev.off()
}


# create plots for targeted genome
for (measure in pathology_measures) {
    pdf(file.path(dir_plots, paste0("spe_targeted", "_", measure, ".pdf")), width = 18)
    for (i in cols_targeted) {
        plot <- ggpubr::ggviolin(
            cluster_targeted_df,
            i,
            measure,
            add = "boxplot",
            add.params = list(fill = "white")
        )
        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE,
            scales = "free_y"
        )
        print(plot)
    }
    dev.off()
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
