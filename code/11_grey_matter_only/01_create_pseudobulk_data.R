library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

genome_type <- c('wholegenome', 'targeted')

k <- 2 #2:28
k_nice <-sprintf("%02d", k)

for(type in genome_type){
    dir.create(here::here("processed-data","gm_only", type), showWarnings = FALSE)

    spe <-
        readRDS(
            here::here(
                "processed-data",
                "08_harmony_BayesSpace",
                type,
                paste0("spe_harmony_",type, ".rds")

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

    # > dim(colData(spe))
    # [1] 38115   108

    spe_new <- spe[,!spe$subject == "Br3874"]
    spe_new <- spe_new[, !spe_new$BayesSpace_harmony_k02 == 2]


    # > dim(colData(spe_new))
    # [1] 21086   108

    ##pseudobulk across pathology labels



}
