# sgejobs::job_loop(
#     loops = list(spefile = c(
#         "spe_postqc.Rdata", "spe_targeted_postqc.Rdata"
#     )),
#     name = "preprocess_and_harmony",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "30G",
#     cores = 4L
# )
# dir.create(here::here("code", "08_harmony_BayesSpace", "logs"), showWarnings = FALSE)

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
if (opt$spefile == "spe_targeted_postqc.Rdata") {
    spe <- spe_targeted
    suffix <- "targeted"
} else {
    suffix <- "wholegenome"
}

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", suffix)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", suffix)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clustering_results"), showWarnings = FALSE)

## Load the data
load(here::here("processed-data", "07_spot_qc", opt$spefile), verbose = TRUE)



message("Running quickCluster()")
set.seed(20220201)
Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(4),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
## Might be needed:
# options(error = recover)
spe <-
    computeSumFactors(spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))

message("Running logNormCounts()")
spe <- logNormCounts(spe)


message("Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(4)
)

pdf(
    file.path(dir_plots, "scran_modelGeneVar.pdf"),
    useDingbats = FALSE
)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message("Running getTopHVGs()")
top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)

save(top.hvgs,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = file.path(dir_rdata, "top.hvgs.Rdata")
)


message("Running runPCA()")
set.seed(20220201)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

pdf(
    file.path(dir_plots, "pca_elbow.pdf"),
    useDingbats = FALSE
)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

message("Running runTSNE() perplexity 5")
Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity05",
        perplexity = 5
    )
Sys.time()

message("Running runTSNE() perplexity 20")
Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity20",
        perplexity = 20
    )
Sys.time()

message("Running runTSNE() perplexity 50")
Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity50",
        perplexity = 50
    )
Sys.time()

message("Running runTSNE() perplexity 80")
Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity80",
        perplexity = 80
    )
Sys.time()


message("Running runUMAP()")
Sys.time()
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()

## Perform harmony batch correction
message("Running RunHarmony()")
Sys.time()
spe <- RunHarmony(spe, "sample_id", verbose = FALSE)
Sys.time()

message("Running runUMAP() on HARMONY dimensions")
Sys.time()
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
Sys.time()
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


## Explore UMAP results
pdf(file = file.path(dir_plots, "UMAP_subject.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))
) +
    geom_point() +
    labs(color = "Subject") +
    theme_bw()
dev.off()

pdf(file = file.path(dir_plots, "UMAP_sample_id.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
    geom_point() +
    labs(color = "sample_id") +
    theme_bw()
dev.off()

## Explore UMAP on HARMONY reduced dimensions
pdf(file = file.path(dir_plots, "UMAP_harmony_sample_id.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
    geom_point() +
    labs(color = "sample_id") +
    theme_bw()
dev.off()


## Perform graph-based clustering on batch corrected-data
message("Running buildSNNGraph() on HARMONY dimensions")
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = "HARMONY")
Sys.time()
save(g_k10, file = file.path(dir_rdata, "g_k10_harmony.Rdata"))

message("Running cluster_walktrap()")
Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
save(g_walk_k10, file = file.path(dir_rdata, "g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
spe$SNN_k10 <- clust_k10 ## Add this one to the SPE too
## Export for later use
cluster_export(
    spe,
    "SNN_k10",
    cluster_dir = file.path(dir_rdata, "clustering_results")
)

message("Running cut_at() from k = 4 to 28")
clust_k5_list <- lapply(4:28, function(n) {
    message(paste(Sys.time(), "n =", n))
    sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0("SNN_k10_k", 4:28)

## Add clusters to spe colData
for (i in seq_along(names(clust_k5_list))) {
    colData(spe) <- cbind(colData(spe), clust_k5_list[i])
    ## Add proper name
    colnames(colData(spe))[ncol(colData(spe))] <- names(clust_k5_list)[i]

    ## Export for later use outside the SPE object
    cluster_export(
        spe,
        names(clust_k5_list)[i],
        cluster_dir = file.path(dir_rdata, "clustering_results")
    )
}

## make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = file.path(dir_plots, "graph_based_harmony.pdf"))
for (i in seq_along(sample_ids)) {
    for (j in seq_along(names(clust_k5_list))) {
        my_plot <- vis_clus(
            spe = spe,
            clustervar = names(clust_k5_list)[j],
            sampleid = sample_ids[i],
            colors = mycolors,
            ... = paste0(" ", names(clust_k5_list)[j])
        )
        print(my_plot)
    }
}
dev.off()


## do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- spatialData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- spatialData(spe)$array_col

pdf(file = file.path(dir_plots, "BayesSpace_offset_check.pdf"))
clusterPlot(spe, "subject", color = NA) + # make sure no overlap between samples
    labs(fill = "Subject", title = "Offset check")
dev.off()


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

## Object size in GB
## (do this near the end in case lobstr crashes, it's happened to me once)
lobstr::obj_size(spe) / 1024^3


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
