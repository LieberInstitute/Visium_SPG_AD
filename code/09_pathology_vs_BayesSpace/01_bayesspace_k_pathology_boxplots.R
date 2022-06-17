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
        "processed-data",
        "08_harmony_BayesSpace",
        "wholegenome",
        "spe_harmony_wholegenome.rds"
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


# create df with relevant variables for whole
cluster_whole_df <- data.frame(
    spot_id = colnames(spe_wholegenome),
    diagnosis = spe_wholegenome$diagnosis,
    sample_id_short = factor(spe_wholegenome$sample_id_short, levels = unique(spe_wholegenome$sample_id_short)),
    NAbeta = spe_wholegenome$NAbeta,
    NpTau = spe_wholegenome$NpTau,
    PAbeta = spe_wholegenome$PAbeta,
    PpTau = spe_wholegenome$PpTau,
    colData(spe_wholegenome)[grep("BayesSpace_harmony_k", colnames(colData(spe_wholegenome)))] ## could use select(matches =) but
    ## have to convert entire colData to df
)


# convert cluster info to factor in cluster_whole_df
cols_whole <- colnames(cluster_whole_df |> select(matches("harmony")))
cluster_whole_df[cols_whole] <- lapply(cluster_whole_df[cols_whole], factor)



# create plots

measure = "PAbeta"
pdf(file.path(dir_plots, paste0("proportion_of_", measure, ".pdf")), width = 18)

for (i in cols_whole) {
    plot <- ggpubr::ggviolin(
        cluster_whole_df,
        i,
        measure,
        add = "boxplot",
        add.params = list(fill = "white"),
        ylab = expression(paste("Proportion of A", beta , " per spot "))
    )
    plot <- facet(plot + theme_bw(),
                  facet.by = "sample_id_short",
                  short.panel.labs = TRUE,
                  scales = "free_y"
    )
    print(plot)
}
dev.off()

measure = "PpTau"

pdf(file.path(dir_plots, paste0("proportion_of_", measure, ".pdf")), width = 18)

for (i in cols_whole) {
    plot <- ggpubr::ggviolin(
        cluster_whole_df,
        i,
        measure,
        add = "boxplot",
        add.params = list(fill = "white"),
        ylab = expression(paste("Proportion of pTau per spot"))
    )
    plot <- facet(plot + theme_bw(),
                  facet.by = "sample_id_short",
                  short.panel.labs = TRUE,
                  scales = "free_y"
    )



    print(plot)
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
