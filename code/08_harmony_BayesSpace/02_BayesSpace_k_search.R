# sgejobs::job_loop(
#     loops = list(spefile = c(
#         "spe_harmony_wholegenome", "spe_harmony_targeted"
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
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")

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
    suffix <- "targeted"
} else if (opt$spefile == "spe_harmony_wholegenome") {
    suffix <- "wholegenome"
}
spe <- readRDS(
    here::here(
        "processed-data", "08_harmony_BayesSpace", suffix,
        paste0("spe_harmony_", suffix, ".rds")
    )
)

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", suffix)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", suffix)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clustering_results"), showWarnings = FALSE)

## Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

set.seed(20220201)

spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 20000)

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_harmony_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clustering_results")
)


spe <- spatialEnhance(spe, use.dimred = "HARMONY", q = k, nrep = 20000, burn.in = 4000)

spe$bayesSpace_enhanced_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_harmony_enhanced_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clustering_results")
)


## Visualize BayesSpace results
sample_ids <- unique(colData(spe)$sample_id)
cols <- Polychrome::palette36.colors(k)
names(cols) <- unique(spe$spatial.cluster)

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_k", k_nice),
    pdf_file = file.path(dir_plots, paste0("BayesSpace_harmony_k", k_nice, ".pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 2
)

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_enhanced_k", k_nice),
    pdf_file = file.path(dir_plots, paste0("BayesSpace_harmony_enhanced_k", k_nice, ".pdf")),
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
