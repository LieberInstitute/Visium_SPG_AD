library("spatialLIBD")
library("lobstr")
library("rafalib")
library("ggplot2")
library("here")
library("sessioninfo")

## Create output directories
dir_rdata <- here("processed-data", "22_gene_pairs_path")
dir_plots <- here("plots", "22_gene_pairs_path")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Read in the data
load(here("code", "05_deploy_app_wholegenome", "spe.Rdata"), verbose = TRUE)
# spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
# lobstr::obj_size(spe)
# 2.29 GB

## Identify pathology group of interest
path_i <- as.integer(Sys.getenv("SGE_TASK_ID"))
if(is.na(path_i)) {
    levels(spe$path_groups)
    # [1] "none"   "Ab"     "n_Ab"   "pTau"   "n_pTau" "both"   "n_both"
    warning("Testing with path_i = 2 which is Ab")
    path_i <- 2
}
path_name <- levels(spe$path_groups)[path_i]

## Subset to AD only for the given pathology group of interest
spe_expr_group <- spe[, spe$diagnosis == "AD" & spe$path_groups == path_name]

## Identify which genes have nonzero counts in each spot
assay(spe_expr_group, "expr") <- counts(spe_expr_group) > 0

## Identify top genes expressed in this pathology group
expr_mean <- rowMeans(assay(spe_expr_group, "expr"))
k <- 300
top_k <- sort(expr_mean, decreasing = TRUE)[seq_len(k)]
head(top_k)
summary(top_k)
## We don't need this object anymore
rm(spe_expr_group)

## Subset to genes with highest mean expressed proportion
spe_expr <- spe[names(top_k), spe$diagnosis == "AD"]
assay(spe_expr, "expr") <- counts(spe_expr) > 0
## Final dimensions
dim(spe_expr)
lobstr::obj_size(spe_expr)

## We no longer need the spe object
rm(spe)

## Compute the combinations of gene pairs
gene_combn <- combn(k, 2)
message("Number of gene combinations: ", ncol(gene_combn))

## Find the combination of genes that are co-expressed
message(Sys.time(), " - computing pair_1")
pair_1 <- assay(spe_expr, "expr")[gene_combn[1, ], ]
lobstr::obj_size(pair_1)
message(Sys.time(), " - computing pair_2")
pair_2 <- assay(spe_expr, "expr")[gene_combn[2, ], ]
lobstr::obj_size(pair_2)
message(Sys.time(), " - computing co_expr")
co_expr <- pair_1 & pair_2
message(Sys.time(), " - done computing co_expr")
lobstr::obj_size(co_expr)

## Build a SummarizedExperiment object with the gene co-expression pairs
co_expr_se <- SummarizedExperiment(
    assays = SimpleList(co_expr = co_expr),
    rowData = DataFrame(
        gene_id_1 = rownames(pair_1),
        gene_id_2 = rownames(pair_2),
        gene_name_1 = rowData(spe_expr)$gene_name[gene_combn[1, ]],
        gene_name_2 = rowData(spe_expr)$gene_name[gene_combn[2, ]]
    ),
    colData = colData(spe_expr)
)
rowData(co_expr_se)$gene_name_combn <- paste0(
    rowData(co_expr_se)$gene_name_1,
    "_",
    rowData(co_expr_se)$gene_name_2
)
## Set the rownames for the gene combination
rownames(co_expr_se) <- paste0(rownames(pair_1), "_", rownames(pair_2))
## Use unique colnames, otherwise the plotting code breaks later on
colnames(co_expr_se) <- co_expr_se$key

message("Memory for co_expr_se:")
lobstr::obj_size(co_expr_se)

## Remove objects we don't need anymore
rm(pair_1, pair_2, co_expr)

## Calculate mean of co-expression across groups of interest
path_list <- rafalib::splitit(spe_expr$path_groups)
co_expr_means <- do.call(cbind, lapply(path_list, function(ii) {
    rowMeans(assay(co_expr_se, "co_expr")[, ii])
}))
dim(co_expr_means)

## Find highest (max) and second highest for each gene pair
rowData(co_expr_se)$highest <- apply(co_expr_means, 1, max)
rowData(co_expr_se)$which_highest <- colnames(co_expr_means)[apply(co_expr_means, 1, which.max)]
rowData(co_expr_se)$second <- apply(co_expr_means, 1, function(x) {
    max(x[x != max(x)])
})
rowData(co_expr_se)$which_second <- apply(co_expr_means, 1, function(x) {
    removed_max <- x[x != max(x)]
    names(which.max(removed_max))
})
rowData(co_expr_se)$ratio <- rowData(co_expr_se)$highest / rowData(co_expr_se)$second
message("Summary of ratio of highest / second highest:")
summary(rowData(co_expr_se)$ratio)

head(sort(rowData(co_expr_se)$ratio, decreasing = TRUE))
# ENSG00000127585_ENSG00000089737 ENSG00000102003_ENSG00000089737 ENSG00000127585_ENSG00000145920
#                        1.285471                        1.278619                        1.266584
# ENSG00000089157_ENSG00000089737 ENSG00000111716_ENSG00000089737 ENSG00000137409_ENSG00000089737
#                        1.258635                        1.257425                        1.251863

message("Full co_expr_se size:")
lobstr::obj_size(co_expr_se)

## Plotting function with code adapted from
## https://github.com/LieberInstitute/Visium_SPG_AD/blob/7f4017344675764da592d27677d34993589e44a5/code/21_spot_deconvo/02_find_markers.R#L71-L122
my_plot_expression <- function(
        sce, genes, assay = "counts", ct = "cellType", title = NULL,
        marker_stats
    ) {
    stopifnot(length(unique(colnames(sce))) == ncol(sce))
    cat_df <- as.data.frame(colData(sce))[, ct, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))
    expression_long$value <- as.integer(expression_long$value)

    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    #   Use gene symbols for labels, not Ensembl ID
    symbols <- rowData(sce)$gene_name_combn[match(genes, rownames(sce))]
    names(symbols) <- genes

    #   Add a data frame for adding mean-ratio labels to each gene
    text_df <- marker_stats[marker_stats$gene %in% genes, ]
    # text_df$ratio <- paste0("Mean ratio: ", round(text_df$ratio, 2))
    text_df$Var1 <- factor(text_df$gene, levels = levels(expression_long$Var1))

    expression_violin <- ggplot(
        data = expression_long, aes(x = cat, y = value, fill = cat)
    ) +
        geom_violin(scale = "width") +
        geom_text(
            data = text_df,
            mapping = aes(
                x = length(unique(sce[[ct]])), y = Inf, fill = NULL,
                label = ratio
            ),
            size = 10, hjust = 1, vjust = 1.3
        ) +
        scale_fill_manual(values = discrete_cell_palette) +
        facet_wrap(
            ~Var1,
            ncol = 5, scales = "free_y",
            labeller = labeller(Var1 = symbols)
        ) +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw(base_size = 35) +
        theme(
            legend.position = "None", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic")
        ) +
        stat_summary(fun = median, geom = "crossbar", width = 0.3) +
        stat_summary(
            fun = "mean",
            geom = "crossbar",
            width = 0.5,
            colour = "red"
        )

    # expression_violin
    return(expression_violin)
}

## Identify the colors for the pathology groups
cols <- unique(co_expr_se$path_groups_colors)
names(cols) <- sapply(cols, function(x) {
    names(co_expr_se$path_groups_colors[co_expr_se$path_groups_colors == x])[1]
})

discrete_cell_palette <- cols


pdf(
    file.path(dir_plots, paste0("gene_pairs_co-expr_", path_name, ".pdf")),
    width = 35, height = 35
)
for(i in seq(0, 90, by = 10)) {
    p <- my_plot_expression(
        sce = co_expr_se,
        genes = names(
            sort(rowData(co_expr_se)$ratio, decreasing = TRUE)
        )[seq_len(10) + i],
        assay = "co_expr",
        ct = "path_groups",
        marker_stats = data.frame(
            ratio = paste0(
                round(rowData(co_expr_se)$ratio, 3),
                " (",
                rowData(co_expr_se)$which_highest,
                "/",
                rowData(co_expr_se)$which_second,
                ")"
            ),
            gene = rownames(co_expr_se)
        )
    )
    print(p)
}
dev.off()

## Save object for later
saveRDS(
    co_expr_se,
    file = file.path(dir_rdata, paste0("co_expr_se_", path_name, ".rds"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
