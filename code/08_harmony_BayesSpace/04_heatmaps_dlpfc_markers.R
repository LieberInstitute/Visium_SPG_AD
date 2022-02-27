library('readr')
library('here')
library('pheatmap')
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("scuttle")
library('stringr')
library(Matrix.utils)
library(ComplexHeatmap)
library(cowplot)




sig_genes <- read_csv('https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/Analysis/Layer_Guesses/sig_genes.csv')
table(sig_genes$test)

## Load basic SPE data
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)


##output directories
dir_plots <- here::here("plots", "07_spot_qc", "heatmaps")
dir.create(dir_plots, showWarnings = FALSE)


#import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") #, suffix

spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")

spe<- logNormCounts(spe)

sig_genes <- sig_genes |> filter(test == "layer_vs_rest")

counts_subset <- counts(spe)[rownames(counts(spe)) %in% sig_genes$ensembl, ]
# counts_subset <-as.data.frame(as.matrix(counts_subset))
sig_genes <- as.data.frame(sig_genes)
rownames(sig_genes) <- sig_genes$ensembl
sig_genes$gene_layer <- paste(sig_genes$gene, sig_genes$layer)
#
# counts_subset <- merge(x = counts_subset, y = sig_genes, by ='row.names')
# counts_subset <- na.omit(counts_subset)
# rownames(counts_subset) <- counts_subset$gene_layer
# counts_subset <- counts_subset[, -which(names(counts_subset) %in% colnames(sig_genes))]
# counts_subset<- select(counts_subset, -1)
# counts_subset <- as(as.matrix(counts_subset), "sparseMatrix")

#assay(spe, "counts_subset") <- counts_subset
#can't add counts_subset to assays since rownames are not identical to OG


rownames(sig_genes)
sig_gen_index <- match(colnames(pb), rownames(sig_genes))
reordered <- sig_genes[sig_gen_index,]


##pseudobulking

spe_new <- SpatialExperiment(assays = list(counts = counts_subset),
                         colData = colData(spe))


groups <- colData(spe_new)[, c("sample_id", "imported_BayesSpace_harmony_k15")]

pb <- aggregate.Matrix(t(counts(spe_new)),
                       groupings = groups, fun = "sum")

pb<- as.matrix(pb)
pb<-as.data.frame(pb)
colnames(pb) <- reordered$gene_layer

pb$sample_cluster <- rownames(pb)
pb <-pb |> mutate(sample = str_sub(sample_cluster, 1,19))
pb<- pb |> mutate(cluster = str_sub(sample_cluster, 21))

s <- split(pb, pb$sample)


## plot for k = 15
myplots <- list()

for(i in 1:10){
rownames(s[[i]]) = s[[i]]$cluster
plot <- grid.grabExpr(draw(Heatmap(as.matrix(s[[i]][,1:69]),
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 6),
        column_title = as.character(s[[i]]$sample[1]),
        column_title_gp = grid::gpar(fontsize = 8),
        name = " ")))
myplots[[i]] <- plot


}
cowplot::plot_grid(plotlist = myplots, ncol =4)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



#ave.beta <- aggregateAcrossCells(spe_new, statistics="mean",
#use.assay.type="counts", ids=c(spe$sample_id, spe$imported_BayesSpace_harmony_k08))
