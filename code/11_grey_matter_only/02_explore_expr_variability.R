# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "02_explore_expr_variability",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 02_explore_expr_variability.sh

# Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}


library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("scater")


## output directory
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE)
dir.create(file.path(dir_rdata, opt$spetype), showWarnings = FALSE)

## load spe data
sce_pseudo <-
    readRDS(
        here::here(
            "processed-data",
            "11_grey_matter_only",
            opt$spetype,
            paste0("sce_pseudo_pathology_", opt$spetype, ".rds")
        )
    )


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
