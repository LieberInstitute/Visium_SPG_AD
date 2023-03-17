library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("sessioninfo")

## Create output directory
dir_rdata <- here("processed-data", "98_prepare_to_share")
dir.create(dir_rdata, showWarnings = FALSE)


## Determine the suffix
suffix <-
    ifelse(as.numeric(Sys.getenv("SGE_TASK_ID")) == 1, "wholegenome", "targeted")

## Create output directory
dir.create(file.path(dir_rdata, suffix), showWarnings = FALSE)

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
            "Ab",
            "n_Ab",
            "pTau",
            "n_pTau",
            "both",
            "n_both"
        )
    )

## Convert pathology variables into factors
for (i in colnames(colData(spe))[grep("^path_", colnames(colData(spe)))]) {
    colData(spe)[[i]] <- factor(colData(spe)[[i]])
}

## Change Braak info based on latest information from LIBD pathology
unique(spe$subject)


unique(spe$BCrating)
unique(spe$braak)
unique(spe$cerad)
spe$BCrating <- NULL ## This variable was removed from the phenotype table
spe$braak <- c("Br3854" = "Stage VI", "Br3873" = "Stage V", "Br3880" = "Stage VI", "Br3874" = "Stage IV")[spe$subject]
spe$cerad <- c("Br3854" = "Frequent", "Br3873" = "Frequent", "Br3880" = "Frequent", "Br3874" = "None")[spe$subject]
unique(spe$braak)
unique(spe$cerad)

## Add APOe genotype info
sce_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E3", "Br3880" = "E3/E3", "Br3874" = "E2/E3")[spe$subject]

## Load pathology colors
## This info is used by spatialLIBD v1.7.18 or newer
source(here("code", "colors_pathology.R"), echo = TRUE, max.deparse.length = 500)
spe$path_groups_colors <- colors_pathology[as.character(spe$path_groups)]

## Save the final object that we can share through spatialLIBD
save(spe, file = file.path(dir_rdata, suffix, "spe.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
