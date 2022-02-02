# sgejobs::job_loop(
#     loops = list(spefile = c(
#         "spe_postqc", "spe_targeted_postqc"
#     )),
#     name = "BayesSpace_k_search",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "80G",
#     task_num = 15
# )

## Required libraries
library("getopt")
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("sessioninfo")
library("scran") ## requires uwot for UMAP
library("uwot")
library("scater")
library("BiocParallel")
library("PCAtools")
library("ggplot2")
library("BayesSpace")

## Specify parameters
spec <- matrix(c(
    "spefile", "s", 2, "character", "SPE file name",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Rename from spe_targeted to spe to simplify the code so it can work with
## either
if (opt$spefile == "spe_harmony_targeted") {
    spe <- readRDS(here::here("processed-data", "08_harmony_BayesSpace", "spe_harmony_targeted.rds"))
    suffix <- "targeted"
} else if (opt$spefile == "spe_harmony_wholegenome") {
    spe <- readRDS(here::here("processed-data", "08_harmony_BayesSpace", "spe_harmony_wholegenome.rds"))
    suffix <- "wholegenome"
}

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", suffix)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", suffix)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clustering_results"), showWarnings = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
