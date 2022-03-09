library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("sessioninfo")

## Determine the suffix
suffix <-
    ifelse(as.numeric(Sys.getenv("SGE_TASK_ID")) == 1, "wholegenome", "targeted")

## Load the data
spe <- readRDS(here::here(
    "processed-data",
    "08_harmony_BayesSpace",
    suffix,
    paste0("spe_harmony_", suffix, ".rds")
))

## Import GraphBased clusters
spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        suffix,
        "clusters_graphbased"
    ),
    prefix = "graph_"
)

## Import BayesSpace clusters
spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        suffix,
        "clusters_BayesSpace"
    ),
    prefix = ""
)

## Import pathology levels
spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
    ),
    prefix = ""
)

## Convert from character to a factor, so they appear in the order
## we want
spe$path_groups <-
    factor(
        spe$path_groups,
        levels = c(
            "none",
            "Ab+",
            "next_Ab+",
            "pT+",
            "next_pT+",
            "both",
            "next_both"
        )
    )

## Convert pathology variables into factors
for (i in colnames(colData(spe))[grep("^path_", colnames(colData(spe)))]) {
    colData(spe)[[i]] <- factor(colData(spe)[[i]])
}


## Save the final object for the shiny app
if (suffix == "wholegenome") {
    save(spe,
        file = here("code", "05_deploy_app_wholegenome", "spe.Rdata")
    )
} else {
    save(spe, file = here("code", "06_deploy_app_targeted", "spe.Rdata"))
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
