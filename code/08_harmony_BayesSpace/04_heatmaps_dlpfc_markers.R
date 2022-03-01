library('readr')
library('here')
library("SpatialExperiment")
library("scran")
library("scater")
library("spatialLIBD")
library("sessioninfo")
library('ComplexHeatmap')
library('cowplot')


##load and process sig_genes
sig_genes <- read_csv('https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/Analysis/Layer_Guesses/sig_genes.csv')
table(sig_genes$test)

#filter layer-related genes
sig_genes <- sig_genes |> filter(test == "layer_vs_rest")
#dim(sig_genes)  70x15

#create new column to merge gene name and layer
sig_genes <- as.data.frame(sig_genes)
rownames(sig_genes) <- sig_genes$ensembl
sig_genes$gene_layer <- paste(sig_genes$gene, sig_genes$layer)

##output directories
dir_plots <- here::here("plots", "07_spot_qc", "heatmaps")
dir.create(dir_plots, showWarnings = FALSE)

################################whole genome####################################

## Load basic SPE data
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)


#import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") #, suffix

spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")

## Check genes not present in our data
sig_genes[!sig_genes$ensembl %in% rownames(spe), ]
#                 top  layer      gene   tstat         pval          fdr gene_index
# ENSG00000259527   4 Layer1 LINC00052 11.3532 4.328627e-18 2.416564e-14      16455
#                         ensembl KM_Zeng    BM RNAscope          test in_rows
# ENSG00000259527 ENSG00000259527   FALSE FALSE     TRUE layer_vs_rest  14;180
#                                         results       gene_layer
# ENSG00000259527 Layer1_top4;Layer1-Layer6_top10 LINC00052 Layer1


#check if BayesSpace columns have NA vals
sum(is.na(as.data.frame(colData(spe)) |> select(matches("BayesSpace_harmony"))))
#[1] 0


bayes_cols <- colnames(as.data.frame(colData(spe)) |> select(matches("BayesSpace_harmony")))


pdf(file.path(dir_plots, paste0("spe_whole_dlpfc_heatmap", ".pdf")), width = 14)
for( k in bayes_cols){
    groups <- colData(spe)[, c("sample_id", k)]

    ## Pseudo-bulk for our current BayesSpace cluster results
    spe_pseudo <- aggregateAcrossCells(
        spe,
        DataFrame(BayesSpace = colData(spe)[[k]], sample_id = spe$sample_id)
    )
    spe_pseudo <- logNormCounts(spe_pseudo)

    ## plot for k = 15
    myplots <- vector("list", 10)

    for(i in 1:10){

        heat_matrix <- logcounts(spe_pseudo)[rownames(spe_pseudo) %in% sig_genes$ensembl, spe_pseudo$sample_id_short == levels(spe_pseudo$sample_id_short)[i] ]
        colnames(heat_matrix) <- unique(spe_pseudo$BayesSpace)
        rownames(heat_matrix) <- sig_genes$gene_layer[match(rownames(heat_matrix), sig_genes$ensembl)]
        rownames(heat_matrix) <- gsub("ayer", "", rownames(heat_matrix))

        plot <- grid.grabExpr(draw(Heatmap(t(heat_matrix),
                                           column_names_gp = grid::gpar(fontsize = 3.2),
                                           row_names_gp = grid::gpar(fontsize = 5),
                                           column_title = spe_pseudo$sample_id[spe_pseudo$sample_id_short == levels(spe_pseudo$sample_id_short)[i]][1],
                                           column_title_gp = grid::gpar(fontsize = 8),
                                           name = " ")))
        myplots[[i]] <- plot


    }
    plot_grid <- cowplot::plot_grid(plotlist = myplots, ncol =4)
    print(plot_grid)
}
dev.off()

################################targeted genome##################################

load(here::here("processed-data", "07_spot_qc", "spe_targeted_postqc.Rdata"), verbose = TRUE)

#import cluster info for targeted
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "targeted") #, suffix

spe_targeted <- cluster_import(
    spe_targeted,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_")

#add logcounts(spe)
spe_targeted<- logNormCounts(spe_targeted)

#only keep gene logcounts found in sig_genes
counts_subset <- logcounts(spe_targeted)[rownames(logcounts(spe_targeted)) %in% sig_genes$ensembl, ]
#dim(counts_subset)
#[1] 68 38115
#setdiff(rownames(counts_subset),sig_genes$ensembl)
#character(0)   setdiff returns this.


##pseudobulking

spe_new <- SpatialExperiment(assays = list(counts = counts_subset),
                             colData = colData(spe_targeted))


bayes_cols <- colnames(as.data.frame(colData(spe_new)) |> select(matches("BayesSpace_harmony")))


pdf(file.path(dir_plots, paste0("spe_targeted_dlpfc_heatmap", ".pdf")), width = 14)
for( k in bayes_cols){
    groups <- colData(spe_new)[, c("sample_id", k)]

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

################################################################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


