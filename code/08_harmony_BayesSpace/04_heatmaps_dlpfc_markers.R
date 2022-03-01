library('readr')
library('here')
library('pheatmap')
library("SpatialExperiment")
library("scran")
library("scater")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library('stringr')
library('Matrix.utils')
library('ComplexHeatmap')
library('cowplot')
library('Polychrome')

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

#add logcounts(spe)
spe<- logNormCounts(spe)

#filter layer-related genes
sig_genes <- sig_genes |> filter(test == "layer_vs_rest")
#dim(sig_genes)  70x15

#create new column to merge gene name and layer
sig_genes <- as.data.frame(sig_genes)
rownames(sig_genes) <- sig_genes$ensembl
sig_genes$gene_layer <- paste(sig_genes$gene, sig_genes$layer)


#only keep gene logcounts found in sig_genes
counts_subset <- logcounts(spe)[rownames(logcounts(spe)) %in% sig_genes$ensembl, ]
#dim(counts_subset)
#[1]    69 38115
#setdiff(rownames(counts_subset),sig_genes$ensembl)
#character(0)   setdiff returns this.





##pseudobulking

spe_new <- SpatialExperiment(assays = list(counts = counts_subset),
                         colData = colData(spe))


#check if BayesSpace columns have NA vals
#sum(is.na(as.data.frame(colData(spe_new)) |> select(matches("BayesSpace_harmony"))))
#returns 0


bayes_cols <- colnames(as.data.frame(colData(spe_new)) |> select(matches("BayesSpace_harmony")))


pdf(file.path(dir_plots, paste0("spe_whole_dlpfc_heatmap", ".pdf")), width = 14)
for( k in bayes_cols){
    groups <- colData(spe)[, c("sample_id", k)]

    pb <- aggregate.Matrix(t(counts(spe_new)),
                       groupings = groups, fun = "sum")

    sig_gen_index <- match(colnames(pb), rownames(sig_genes))
    reordered <- sig_genes[sig_gen_index,]


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
                                           column_names_gp = grid::gpar(fontsize = 3.2),
                                           row_names_gp = grid::gpar(fontsize = 5),
                                           column_title = as.character(s[[i]]$sample[1]),
                                           column_title_gp = grid::gpar(fontsize = 8),
                                           name = " ")))
        myplots[[i]] <- plot


    }
    plot_grid <- cowplot::plot_grid(plotlist = myplots, ncol =4)
    print(plot_grid)
}
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


