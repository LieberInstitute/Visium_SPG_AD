##import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("tidyr")
library("ggpubr")
library("ggplot2")


## Load basic SPE data
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)


##output directories
dir_plots <- here::here("plots", "07_spot_qc", "pathology_vs_Bayesspace_cluster_barplots")
dir.create(dir_plots, showWarnings = FALSE)


#import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") #, suffix

cluster_spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")



controls <- c("V10A27106_A1_Br3874", "V10A27004_A1_Br3874", "V10T31036_A1_Br3874")
#last one shouldn't be used for pTau


colData(cluster_spe)
#imported_BayesSpace_harmony_k15

cluster_df <- as.data.frame(colData(cluster_spe))

cluster_df <- cluster_df |> mutate(across(matches("BayesSpace_harmony"),factor ))

cluster_df <- cluster_df |> mutate(pTau_outliers =
                                       ifelse(NpTau > 7 &PpTau > 0.014, "outlier", "normal"))



# Stacked bar plots, add labels inside bars
plot <- ggbarplot(cluster_df, x = "imported_BayesSpace_harmony_k14", y = "pTau_outliers",
          fill = "pTau_outliers", color = "pTau_outliers",
          palette = c("gray", "black"),
          label = TRUE, lab.col = "white", lab.pos = "in")

