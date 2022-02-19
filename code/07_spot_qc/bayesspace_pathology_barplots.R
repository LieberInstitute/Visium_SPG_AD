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






colData(cluster_spe)
#imported_BayesSpace_harmony_k15

cluster_df <- as.data.frame(colData(cluster_spe))

cluster_df <- cluster_df |> mutate(across(matches("BayesSpace_harmony"),factor ))

cluster_df <- cluster_df |> mutate(pTau_outliers =
                                       ifelse(NpTau > 7 & PpTau > 0.014, "outlier", "normal"))



# Stacked bar plots, add labels inside bars

plot <- ggplot(cluster_df, aes(fill = pTau_outliers,
                       x = imported_BayesSpace_harmony_k14, y = PpTau))+
    geom_bar(position = "fill", stat = "identity")



bayes_cols <- cluster_df |> select(matches("BayesSpace_harmony"))

##

pdf(file.path(dir_plots, paste0("spe_whole_barplots_pTau", ".pdf")), width = 14)

for (i in bayes_cols) {
    plot <- ggplot(cluster_df, aes(x= i))+
        geom_bar(aes(fill = pTau_outliers),
                    position = 'fill')+
        scale_y_continuous(labels = scales::percent)+
        labs ( y = "Percentage")
    plot <-facet(plot+ theme_bw(), facet.by = "sample_id",
            short.panel.labs = FALSE)


    print(plot)
    }
dev.off()



