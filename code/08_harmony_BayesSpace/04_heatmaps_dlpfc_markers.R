library("readr")
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("spatialLIBD")
library("sessioninfo")
library("ComplexHeatmap")
library("cowplot")
library("dplyr")


## load and process sig_genes
sig_genes <-
    read_csv(
        "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/Analysis/Layer_Guesses/sig_genes.csv"
    )
table(sig_genes$test)

# filter layer-related genes
sig_genes <- sig_genes |> dplyr::filter(test == "layer_vs_rest")
# dim(sig_genes)  70x15

# create new column to merge gene name and layer
sig_genes <- as.data.frame(sig_genes)
rownames(sig_genes) <- sig_genes$ensembl
sig_genes$gene_layer <- paste(sig_genes$gene, sig_genes$layer)

## output directories
dir_plots <-
    here::here("plots", "08_harmony_BayesSpace", "marker_heatmaps")
dir.create(dir_plots, showWarnings = FALSE)


## Function for making the heatmap

heatmap_spe <- function(suffix) {
    ## Load basic SPE data
    spe <- readRDS(here::here(
        "processed-data",
        "07_spot_qc",
        paste0("spe_", suffix, "_postqc.rds")
    ))

    # import cluster info
    spe <- cluster_import(
        spe,
        cluster_dir = here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            suffix,
            "clusters_BayesSpace"
        ),
        prefix = ""
    )

    ## Check genes not present in our data
    print(sig_genes[!sig_genes$ensembl %in% rownames(spe), ])
    #                 top  layer      gene   tstat         pval          fdr gene_index
    # ENSG00000259527   4 Layer1 LINC00052 11.3532 4.328627e-18 2.416564e-14      16455
    #                         ensembl KM_Zeng    BM RNAscope          test in_rows
    # ENSG00000259527 ENSG00000259527   FALSE FALSE     TRUE layer_vs_rest  14;180
    #                                         results       gene_layer
    # ENSG00000259527 Layer1_top4;Layer1-Layer6_top10 LINC00052 Layer1


    # check if BayesSpace columns have NA vals
    print(sum(is.na(
        as.data.frame(colData(spe)) |> select(matches("BayesSpace_harmony"))
    )))
    # [1] 0


    bayes_cols <-
        colnames(as.data.frame(colData(spe)) |> select(matches("BayesSpace_harmony")))


    pdf(file.path(
        dir_plots,
        paste0("spe_", suffix, "_dlpfc_heatmap", ".pdf")
    ), width = 14)
    for (k in bayes_cols) {
        groups <- colData(spe)[, c("sample_id", k)]

        ## Pseudo-bulk for our current BayesSpace cluster results
        spe_pseudo <- aggregateAcrossCells(
            spe,
            DataFrame(
                BayesSpace = colData(spe)[[k]],
                sample_id = spe$sample_id
            )
        )
        spe_pseudo <- logNormCounts(spe_pseudo)

        ## plot for k = 15
        myplots <- vector("list", 10)

        for (i in 1:10) {
            heat_matrix <-
                logcounts(spe_pseudo)[rownames(spe_pseudo) %in% sig_genes$ensembl, spe_pseudo$sample_id_short == levels(spe_pseudo$sample_id_short)[i]]
            colnames(heat_matrix) <- unique(spe_pseudo$BayesSpace)
            rownames(heat_matrix) <-
                sig_genes$gene_layer[match(rownames(heat_matrix), sig_genes$ensembl)]
            rownames(heat_matrix) <-
                gsub("ayer", "", rownames(heat_matrix))

            plot <- grid.grabExpr(draw(
                Heatmap(
                    t(heat_matrix),
                    column_names_gp = grid::gpar(fontsize = 3.2),
                    row_names_gp = grid::gpar(fontsize = 5),
                    column_title = spe_pseudo$sample_id[spe_pseudo$sample_id_short == levels(spe_pseudo$sample_id_short)[i]][1],
                    column_title_gp = grid::gpar(fontsize = 8),
                    name = " "
                )
            ))
            myplots[[i]] <- plot
        }
        plot_grid <-
            cowplot::plot_grid(plotlist = myplots, ncol = 4)
        print(plot_grid)
    }
    dev.off()
}

## Run the function
heatmap_spe("wholegenome")
heatmap_spe("targeted")

################################################################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
