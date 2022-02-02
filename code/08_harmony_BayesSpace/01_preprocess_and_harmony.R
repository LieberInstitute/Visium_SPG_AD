## Required libraries
library("getopt")
library("here")
library("SpatialExperiment")
library("sessioninfo")
library("scran") ## requires uwot for UMAP
library("uwot")
library("scater")
library("BiocParallel")
library("PCAtools")

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


set.seed(20220201)
Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(4),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()

Sys.time()
# [1] "2021-12-15 13:19:48 EST"
## Might be needed:
# options(error = recover)
spe <-
    computeSumFactors(spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

table(spe$scran_quick_cluster)

summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000034  0.477444  0.847973  1.000000  1.331274 16.813313

spe <- logNormCounts(spe)


## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(4)
)

pdf(
    file.path(dir_plots, "scran_modelGeneVar_final.pdf"),
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

top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 18417

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 17830

save(top.hvgs,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = file.path(dir_rdata, "top.hvgs_all_final.Rdata")
)



set.seed(20220201)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
# 50

pdf(
    file.path(dir_plots, "pca_elbow_final.pdf"),
    useDingbats = FALSE
)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity50",
        perplexity = 50
    )
Sys.time()
# [1] "2021-02-17 10:45:30 EST"
# [1] "2021-02-17 11:02:59 EST"


Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity5",
        perplexity = 5
    )
Sys.time()
# [1] "2021-09-08 11:14:41 EDT"
# [1] "2021-09-08 12:03:04 EDT"

Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity20",
        perplexity = 20
    )
Sys.time()
# [1] "2021-09-08 12:17:38 EDT"
# [1] "2021-09-08 13:09:17 EDT"

Sys.time()
set.seed(20220201)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity80",
        perplexity = 80
    )
Sys.time()




# RunUMAP
Sys.time()
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()





# make plots ofUMAP
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

### harmony batch correction
spe <- RunHarmony(spe, "sample_id", verbose = F)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


pdf(file = file.path(dir_plots, "UMAP_harmony_sample_id.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
    geom_point() +
    labs(color = "sample_id") +
    theme_bw()
dev.off()


## graph-based on batch correct
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = "HARMONY")
Sys.time()
# "2021-12-16 16:17:31 EST"
# "2021-12-16 16:34:28 EST"
save(g_k10, file = file.path(dir_rdata, "g_k10_harmony.Rdata"))

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
# [1] "2021-12-17 13:05:59 EST"

save(g_walk_k10, file = file.path(dir_rdata, "g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = file.path(dir_rdata, "clust_k10_harmony.Rdata"))

clust_k5_list <- lapply(4:28, function(n) {
    message(paste(Sys.time(), "n =", n))
    sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0("k", 4:28)
save(clust_k5_list, file = file.path(dir_rdata, "clust_k5_list_harmony.Rdata"))

## Add clusters to spe colData
cluster_colNames <- paste0("SNN_k10_k", 4:28)
for (i in seq_along(cluster_colNames)) {
    colData(spe) <- cbind(colData(spe), clust_k5_list[i])
}
colnames(colData(spe))[34:58] <- cluster_colNames

## make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = file.path(dir_plots, "vis_clus_graph_based_har.pdf"))
for (i in seq_along(sample_ids)) {
    for (j in seq_along(cluster_colNames)) {
        my_plot <- vis_clus(
            spe = spe,
            clustervar = cluster_colNames[j],
            sampleid = sample_ids[i],
            colors = mycolors,
            ... = paste0(" ", cluster_colNames[j])
        )
        print(my_plot)
    }
}
dev.off()

cluster_export(
    spe,
    "SNN_k10_k7",
    cluster_dir = file.path(dir_rdata, "clustering_results")
)

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


## do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- spatialData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- spatialData(spe)$array_col

pdf(file = file.path(dir_plots, "bayesSpace_offset_check.pdf"))
clusterPlot(spe, "subject", color = NA) + # make sure no overlap between samples
    labs(fill = "Subject", title = "Offset check")
dev.off()


## Object size in GB
lobstr::obj_size(spe) / 1024^3


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
