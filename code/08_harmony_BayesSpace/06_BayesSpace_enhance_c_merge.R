# sgejobs::job_single(
#     name = "BayesSpace_enhance_c_merge",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "30G",
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
library("SingleCellExperiment")
library("BayesSpace")
library("Polychrome")


## create output directories
dir_rdata <-
    here::here("processed-data", "08_harmony_BayesSpace", opt$spetype)
dir_plots <- here::here("plots", "08_harmony_BayesSpace", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Locate RDS files
rds_files <-
    dir(
        dir_rdata,
        pattern = paste0("^sce_enhanced_", opt$spetype, "_sample"),
        full.names = TRUE
    )

## Load the data and combine into a single object
sce_enhanced <- do.call(cbind, lapply(rds_files, readRDS))

## Manually pass what k value we selected with fasthplus
k_selected <- ifelse(opt$spetype == "wholegenome", 2, 2)

## Make a plot of the enhanced cluster results
pdf(file.path(dir_plots,
    "BayesSpace_enhanced.pdf"
), width = 7 * length(rds_files))
clusterPlot(sce_enhanced,
    palette = setNames(
        Polychrome::palette36.colors(k_selected),
        seq_len(k_selected)
    ))
dev.off()

## Save for later
saveRDS(sce_enhanced, file = file.path(dir_rdata, file.path(
    dir_rdata,
    paste0("sce_enhanced_",
        opt$spetype,
        ".rds")
)))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
