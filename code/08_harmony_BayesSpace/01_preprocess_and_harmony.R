## Required libraries
library("getopt")
library("here")
library("SpatialExperiment")
library("sessioninfo")

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

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace")
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load the data
load(here::here("processed-data", "07_spot_qc", opt$spefile), verbose = TRUE)

## Rename from spe_targeted to spe to simplify the code so it can work with
## either
if (opt$spefile == "spe_targeted_postqc.Rdata") {
    spe <- spe_targeted
}





## Save new SPE objects
if (opt$spefile == "spe_targeted_postqc.Rdata") {
    spe_targeted <- spe
    ## First time switching the order of the keywords: now targeted is at the
    ## end, which will make it easier to sort the spe files.
    saveRDS(spe_targeted, file = file.path(dir_rdata, "spe_harmony_targeted.rds"))
} else {
    ## First time using "wholegenome" in the spe name, to clearly differentiate
    ## it from the "targeted" one
    spe_wholegenome <- spe
    saveRDS(spe_wholegenome, file.path(dir_rdata, "spe_harmony_wholegenome.rds"))
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
