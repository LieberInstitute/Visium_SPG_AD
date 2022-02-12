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

##output directories
dir_plots <- here::here("plots", "07_spot_qc", "temp")
dir.create(dir_plots, showWarnings = FALSE)


#import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") #, suffix

cluster_spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")

#import cluster info for targeted genome
dir_rdata_targeted <- here::here("processed-data", "08_harmony_BayesSpace", "targeted")

cluster_spe_targeted <- cluster_import(
    spe_targeted,
    cluster_dir = file.path(dir_rdata_targeted, "clusters_BayesSpace"),
    prefix = "imported_")


#create df with relevant variables for whole
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

#create df with relevant variables for targeted
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


#convert cluster info to factor in cluster_whole_df
cols_whole <- colnames(cluster_whole_df |> select(matches("harmony")))
cluster_whole_df[cols_whole] <- lapply(cluster_whole_df[cols_whole], factor)

#convert cluster info to factor in cluster_targeted_df
cols_targeted <- colnames(cluster_targeted_df |> select(matches("harmony")))
cluster_targeted_df[cols_targeted] <- lapply(cluster_targeted_df[cols_targeted], factor)



pathology_measures = c("PpTau", "PAbeta")

#create plots for whole genome
for( measure in pathology_measures){
    for (i in cols_whole) {
        pdf(file.path(dir_plots, paste0(spe_whole,"_", sample_id,"_", i, ".pdf")), width = 14)

        plot<- ggpubr::ggviolin(
            cluster_whole_df,
            i,
            "PpTau",
            add= "boxplot",
            add.params = list(fill = "white"))
        print(plot)
    }

}



#create plots for targeted genome
for( measure in pathology_measures){
    for (i in cols_targeted) {
        pdf(file.path(dir_plots, paste0(spe_targeted,"_", sample_id,"_", i, ".pdf")), width = 14)

        plot<- ggpubr::ggviolin(
            cluster_targeted_df,
            i,
            "PpTau",
            add= "boxplot",
            add.params = list(fill = "white"))
        print(plot)
    }

}



