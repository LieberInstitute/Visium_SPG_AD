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
library(Matrix.utils)

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

sig_genes <- sig_genes |> filter(test == "layer_vs_rest")
counts_subset <- counts(spe)[rownames(counts(spe)) %in% sig_genes$ensembl, ]

#assay(spe, "counts_subset") <- counts_subset
#can't add counts_subset to assays since rownames are not identical to OG

##pseudobulking


spe_new <- SummarizedExperiment(assays = list(counts = counts_subset),
                         colData = colData(spe))
                         #rowData = rowData(spe)[colnames(rowData(spe)) %in% sig_genes$ensembl] )

groups <- colData(spe_new)[, c("sample_id", "imported_BayesSpace_harmony_k15")]

pb <- aggregate.Matrix(t(counts(spe_new)),
                       groupings = groups, fun = "sum")

#ave.beta <- aggregateAcrossCells(spe_new, statistics="mean",
                                 #use.assay.type="counts", ids=c(spe$sample_id, spe$imported_BayesSpace_harmony_k08))


