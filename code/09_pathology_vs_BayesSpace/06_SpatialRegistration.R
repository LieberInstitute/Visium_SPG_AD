library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)

k <- 4 #as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

#load post BayesSpace spe object
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            "wholegenome",
            "spe_harmony_wholegenome.rds"
        )
    )

spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data",
                             "08_harmony_BayesSpace",
                             "wholegenome",
                             "clusters_BayesSpace"),
    prefix = ""
)


## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        BayesSpace = colData(spe)[[k]],
        sample_id = spe$sample_id
    )
)

#spe_pseudo <- logNormCounts(spe_pseudo) #size factors <0?

# range(spe$size_factor)
# [1]  Inf -Inf
# spe$size_factor[spe$size_factor =='NA'] #size factors
# NULL
#


mat <- counts(spe_pseudo)

# remove lowly expressed genes?
# gIndex = rowMeans(mat) > 0.2
# mat_filter = mat[gIndex, ]


