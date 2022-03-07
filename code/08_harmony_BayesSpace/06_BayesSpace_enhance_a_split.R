# sgejobs::job_single(
#     name = "BayesSpace_enhance_a_split",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "10G",
#     task_num = 2
# )

## Use a simple array job in this case
if (as.numeric(Sys.getenv("SGE_TASK_ID")) == 1) {
    opt <- list(spetype = "wholegenome")
} else {
    opt <- list(spetype = "targeted")
}

## Load required packages
library("here")
library("sessioninfo")
library("SpatialExperiment")

## create output directories
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load the data
spe <- readRDS(here::here(
    "processed-data",
    "08_harmony_BayesSpace",
    opt$spetype,
    paste0("spe_harmony_", opt$spetype, ".rds")
))

## Import BayesSpace cluster initial results
spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        opt$spetype,
        "clusters_BayesSpace"
    ),
    prefix = ""
)

## Manually pass what k value we selected with fasthplus
k_nice_selected <- ifelse(opt$spetype == "wholegenome", "02", "02")

## Rename for BayesSpace
spe$spatial.cluster <-
    colData(spe)[[paste0("BayesSpace_harmony_k", k_nice_selected)]]

## Extra information required by BayesSpace
spe$imagerow <- spe$array_row
spe$imagecol <- spe$array_col

## Remove the image data since it's not used by BayesSpace. This should help
## lower the memory usage
imgData(spe) <- NULL

## Save small spe objects: one for each sample
for (i in seq_len(length(unique(spe$sample_id)))) {
    sample <- unique(spe$sample_id)[i]
    message(Sys.time(), " processing sample ", sample)
    saveRSE(spe[, spe$sample_id == sample], file = file.path(
        dir_rdata,
        paste0("spe_", opt$spetype, "_sample", sprintf("%02d", i), ".rds")
    ))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
