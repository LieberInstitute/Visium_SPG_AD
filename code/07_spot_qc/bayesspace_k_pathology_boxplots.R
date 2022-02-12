##import required libraries
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


dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") #, suffix

cluster_spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")

dir_rdata_targeted <- here::here("processed-data", "08_harmony_BayesSpace", "targeted")

cluster_spe_targeted <- cluster_import(
    spe_targeted,
    cluster_dir = file.path(dir_rdata_targeted, "clusters_BayesSpace"),
    prefix = "imported_")

cluster_whole_df <- data.frame(
    spot_id =rownames(colData(cluster_spe)),
    diagnosis = colData(cluster_spe)$diagnosis,
    sample_id = colData(cluster_spe)$sample_id,
    NAbeta = colData(cluster_spe)$NAbeta,
    NpTau = colData(cluster_spe)$NpTau,
    PAbeta = colData(cluster_spe)$PAbeta,
    PpTau = colData(cluster_spe)$PpTau,
    colData(cluster_spe)[49:59] ##could use select(matches =) but
                                ## have to convert entire colData to df
)

cluster_targeted_df <- data.frame(
    spot_id =rownames(colData(cluster_spe_targeted)),
    diagnosis = colData(cluster_spe_targeted)$diagnosis,
    sample_id = colData(cluster_spe_targeted)$sample_id,
    NAbeta = colData(cluster_spe_targeted)$NAbeta,
    NpTau = colData(cluster_spe_targeted)$NpTau,
    PAbeta = colData(cluster_spe_targeted)$PAbeta,
    PpTau = colData(cluster_spe_targeted)$PpTau,
    colData(cluster_spe_targeted)[49:60]
)


cols <- colnames(cluster_whole_df |> select(matches("harmony")))

cluster_whole_df[cols] <- lapply(cluster_whole_df[cols], factor)

plot_list <- list()

for (i in (cluster_whole_df)|> select(matches("harmony")) ){
    plot<- ggplot(cluster_whole_df, aes( x = i, y = PpTau)) +geom_boxplot()+
    xlab("Cluster ") + ylab("Percentage of pTau")
    plot_list[[colnames(i)]] <- plot
    }






