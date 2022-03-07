# sgejobs::job_loop(
#     loops = list(spetype = c(
#         "wholegenome", "targeted"
#     )),
#     name = "BayesSpace_k_search",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "80G",
#     task_num = 15
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

## Load remaining required packages
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")


## Load the data
spe <- readRDS(
    here::here(
        "processed-data", "08_harmony_BayesSpace", opt$spetype,
        paste0("spe_harmony_", opt$spetype, ".rds")
    )
)

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", opt$spetype)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", opt$spetype)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clusters_BayesSpace"), showWarnings = FALSE)

## Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(20220201)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k)
Sys.time()

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_harmony_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace")
)

## Visualize BayesSpace results
sample_ids <- unique(spe$sample_id)
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_k", k_nice),
    pdf_file = file.path(dir_plots, paste0("BayesSpace_harmony_k", k_nice, ".pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 2
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
