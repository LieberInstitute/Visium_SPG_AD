# sgejobs::job_loop(
#     loops = list(spetype = c(
#         "wholegenome", "targeted"
#     )),
#     name = "BayesSpace_enhance_b_compute",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "50G",
#     task_num = 10
# )

## Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "spetype", "s", 2, "character", "SPE type: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Load required packages
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("BayesSpace")

## create output directories
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Choose what sample we are working with
i <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## Load the data
spe <- readRDS(file.path(
    dir_rdata,
    paste0("spe_", opt$spetype, "_sample", sprintf("%02d", i), ".rds")
))

## Manually pass what k value we selected with fasthplus
k_selected <- ifelse(opt$spetype == "wholegenome", 2, 2)

## Run spatialEnhance for the chosen sample
set.seed(20220307)
sce_enhanced <- spatialEnhance(spe, use.dimred = "HARMONY", q = k_selected)

## Save the results
saveRDS(sce_enhanced, file = file.path(dir_rdata, file.path(
    dir_rdata,
    paste0(
        "sce_enhanced_",
        opt$spetype,
        "_sample",
        sprintf("%02d", i),
        ".rds"
    )
)))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
