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

## Code crashed during spatialEnhance()
## Load the current results and make plots
## then try to debug spatialEnhance()
# spe <- cluster_import(spe, cluster_dir = file.path(dir_rdata, "clusters_BayesSpace"), prefix = "")
# for(k in 4:14) {
#     k_nice <- sprintf("%02d", k)
#     sample_ids <- unique(spe$sample_id)
#     cols <- Polychrome::palette36.colors(k)
#     names(cols) <- sort(unique(colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]]))
#
#     vis_grid_clus(
#         spe = spe,
#         clustervar = paste0("BayesSpace_harmony_k", k_nice),
#         pdf_file = file.path(dir_plots, paste0("BayesSpace_harmony_k", k_nice, ".pdf")),
#         sort_clust = FALSE,
#         colors = cols,
#         spatial = FALSE,
#         point_size = 2
#     )
# }
## For debugging spatialEnhance()
# k <- 4
# k_nice <- "04"
# spe$spatial.cluster <- colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]]

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


message("Running spatialEnhance() -- currently crashes due to https://github.com/edward130603/BayesSpace/issues/71")
Sys.time()
set.seed(20220203)
spe$imagerow <- spatialData(spe)$array_row
spe$imagecol <- spatialData(spe)$array_col

for(sample in unique(spe$sample_id)) {
    spe_small <- spatialEnhance(spe[, spe$sample_id == sample], use.dimred = "HARMONY", q = k)
    spe$spatial.cluster[spe$sample_id == sample] <- spe_small$spatial.cluster
    rm(spe_small)
}
Sys.time()

spe$bayesSpace_enhanced_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_harmony_enhanced_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace")
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
