library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("sessioninfo")

## Determine the suffix
suffix <-
    ifelse(as.numeric(Sys.getenv("SGE_TASK_ID")) == 1, "wholegenome", "targeted")

## Load the data
load(here::here(
    "processed-data",
    "98_prepare_to_share",
    paste0("Visium_SPG_AD_spe_", suffix, ".Rdata")
), verbose = TRUE)

## Drop images we don't really use in the app
imgData(spe) <- imgData(spe)[
    !imgData(spe)$image_id %in% c("hires", "detected", "aligned"),
]

## Add in the spot deconvolution results
spe <- cluster_import(spe, here('processed-data', '21_spot_deconvo'), prefix = 'c2l_')
spe$c2l_sample <- NULL

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
