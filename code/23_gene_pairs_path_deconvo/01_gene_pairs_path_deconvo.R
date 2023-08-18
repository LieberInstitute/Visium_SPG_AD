library("spatialLIBD")
library("lobstr")
library("rafalib")
library("ggplot2")
library("colorspace")
library("here")
library("sessioninfo")

## Create output directories
dir_rdata <- here("processed-data", "23_gene_pairs_path_deconvo")
dir_plots <- here("plots", "23_gene_pairs_path_deconvo")
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

## Subset to just gray matter spots only in the AD samples
spe <- spe[, spe$BayesSpace_harmony_k02 != 2 & spe$diagnosis == "AD"]
# lobstr::obj_size(spe)
# 756.83 MB

## Add in the spot deconvolution results
spe <- cluster_import(spe, here('processed-data', '21_spot_deconvo'), prefix = 'c2l_')
spe$c2l_sample <- NULL

## Get the cell2location column names
c2l_vars <- colnames(colData(spe))[grep("^c2l_", colnames(colData(spe)))]
c2l_vars
# [1] "c2l_ast" "c2l_ex"  "c2l_in"  "c2l_mic" "c2l_oli"
# [6] "c2l_opc"
summary(as.data.frame(colData(spe)[, c2l_vars]))
#    c2l_ast            c2l_ex            c2l_in
# Min.   :0.02672   Min.   :0.00895   Min.   :0.01534
# 1st Qu.:0.14204   1st Qu.:0.13490   1st Qu.:0.11869
# Median :0.21631   Median :0.24416   Median :0.18141
# Mean   :0.25245   Mean   :0.29233   Mean   :0.21439
# 3rd Qu.:0.32461   3rd Qu.:0.39566   3rd Qu.:0.27351
# Max.   :1.48163   Max.   :2.02673   Max.   :1.47782
#    c2l_mic           c2l_oli            c2l_opc
# Min.   :0.03808   Min.   :0.004683   Min.   :0.02605
# 1st Qu.:0.16944   1st Qu.:0.051622   1st Qu.:0.12517
# Median :0.21673   Median :0.095126   Median :0.17422
# Mean   :0.23364   Mean   :0.126026   Mean   :0.20105
# 3rd Qu.:0.27773   3rd Qu.:0.167192   3rd Qu.:0.24741
# Max.   :1.14906   Max.   :0.961844   Max.   :1.10403
c2l_combn <- combn(c2l_vars, 2)
colnames(c2l_combn) <- apply(c2l_combn, 2, function(x) paste0(gsub("c2l_", "", x), collapse = "_"))

cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 1 & colData(spe)[[x[2]]] >= 1
})
summary(rowSums(cell_combn))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.000000 0.000000 0.000332 0.000000 1.000000
sort(table(rowSums(cell_combn)))
# 1     0
# 7 21079
## Almost no spots have 2 cells from cell2location with counts >= 1

summary(rowMeans(as.data.frame(colData(spe)[, c2l_vars])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02297 0.15346 0.20728 0.21998 0.27175 0.89151

## Try with 0.5 as the threshold
cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 0.5 & colData(spe)[[x[2]]] >= 0.5
})
sort(table(rowSums(cell_combn)))
# 15    10     6     3     1     0
#  5    20    48   171   782 20060
table(rowSums(cell_combn[rowSums(cell_combn) > 0, ]) == 1)
# FALSE  TRUE
#   244   782
## 782 spots are specific combinations of just 2 cells (244 have more than 2)

## Try with 0.25 as the threshold
cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 0.25 & colData(spe)[[x[2]]] >= 0.25
})
sort(table(rowSums(cell_combn)))
#  15    10     6     3     1     0
# 447  1475  1971  2732  4139 10322
table(rowSums(cell_combn[rowSums(cell_combn) > 0, ]) == 1)
# FALSE  TRUE
#  6625  4139
## Ok, there's more than with 0.5 that are specific

## Try with 0.15 as the threshold
cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 0.15 & colData(spe)[[x[2]]] >= 0.15
})
sort(table(rowSums(cell_combn)))
#    1    0    3   15    6   10
# 2357 2642 3107 3443 3706 5831
table(rowSums(cell_combn[rowSums(cell_combn) > 0, ]) == 1)
# FALSE  TRUE
# 16087  2357
## It got much noisier than with 0.25

## Try with 0.2 as the threshold
cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 0.2 & colData(spe)[[x[2]]] >= 0.2
})
sort(table(rowSums(cell_combn)))
#   15   10    6    3    1    0
# 1283 3171 3280 3369 3615 6368
table(rowSums(cell_combn[rowSums(cell_combn) > 0, ]) == 1)
# FALSE  TRUE
# 11103  3615
## It' got much noisier than's still noisier than with 0.25

## Let's use 0.25 for now then
cell_combn <- apply(c2l_combn, 2, function(x) {
    colData(spe)[[x[1]]] >= 0.25 & colData(spe)[[x[2]]] >= 0.25
})
spe$c2l_combn_is_unique <- rowSums(cell_combn) == 1
spe$c2l_combn_unique <- rep(NA, ncol(spe))
spe$c2l_combn_unique[spe$c2l_combn_is_unique] <- colnames(cell_combn)[apply(cell_combn[spe$c2l_combn_is_unique, ], 1, which.max)]
sort(table(spe$c2l_combn_unique, useNA = "ifany"))
# oli_opc  in_oli mic_opc mic_oli  ex_oli  in_opc
#      33      40      86      95     106     106
# ast_oli  in_mic ast_opc  ast_in  ex_opc ast_mic
#     153     158     257     289     311     400
#  ex_mic   ex_in  ast_ex    <NA>
#     584     633     888   16947
table(spe$c2l_combn_unique, spe$path_groups, useNA = "ifany")
#         none   Ab n_Ab pTau n_pTau both n_both
# ast_ex   125   47  127  387    116   36     50
# ast_in    42   18   27  146     36    8     12
# ast_mic   73   22   44  160     76   12     13
# ast_oli   44    1    4   45     48    2      9
# ast_opc   40    5   30  100     63    5     14
# ex_in     73   48   78  301     81   18     34
# ex_mic   101   31   64  261     74   22     31
# ex_oli    32    3    6   31     26    2      6
# ex_opc    34   16   32  156     39   14     20
# in_mic    26    4    9   90     20    6      3
# in_oli     7    1    2   20      4    3      3
# in_opc     8    5    9   58     14    3      9
# mic_oli   24    5    4   32     25    0      5
# mic_opc   13    2    6   39     18    5      3
# oli_opc    5    5    0    9     12    0      2
# <NA>    3227  859 1842 6975   2565  591    888

table(spe$c2l_combn_unique, spe$path_groups, useNA = "ifany") >= 10
#         none    Ab  n_Ab  pTau n_pTau  both n_both
# ast_ex   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE
# ast_in   TRUE  TRUE  TRUE  TRUE   TRUE FALSE   TRUE
# ast_mic  TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE
# ast_oli  TRUE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# ast_opc  TRUE FALSE  TRUE  TRUE   TRUE FALSE   TRUE
# ex_in    TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE
# ex_mic   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE
# ex_oli   TRUE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# ex_opc   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE
# in_mic   TRUE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# in_oli  FALSE FALSE FALSE  TRUE  FALSE FALSE  FALSE
# in_opc  FALSE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# mic_oli  TRUE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# mic_opc  TRUE FALSE FALSE  TRUE   TRUE FALSE  FALSE
# oli_opc FALSE FALSE FALSE FALSE   TRUE FALSE  FALSE
# <NA>     TRUE  TRUE  TRUE  TRUE   TRUE  TRUE   TRUE

## Subset to the given pathology group of interest
## + having a unique combination of 2 cells
## + having at least 10 spots where this is the case
spe <- spe[, spe$path_groups == path_name &
            spe$c2l_combn_is_unique]
## + having at least 10 spots where this is the case
combn_has_min10 <- names(table(spe$c2l_combn_unique)[table(spe$c2l_combn_unique) >= 10])
spe <- spe[, spe$c2l_combn_unique %in% combn_has_min10]

## Identify which genes have nonzero counts in each spot
assay(spe, "expr") <- counts(spe) > 0
lobstr::obj_size(spe)
# 197.66 MB

## Identify top genes expressed in this pathology group
expr_mean <- rowMeans(assay(spe, "expr"))
k <- 1000
top_k <- sort(expr_mean, decreasing = TRUE)[seq_len(k)]
head(top_k)
summary(top_k)

## Subset to genes with highest mean expressed proportion
spe_expr <- spe[names(top_k), ]
## Final dimensions
dim(spe_expr)
# [1]   1000 213

## Compute the combinations of gene pairs
gene_combn <- combn(k, 2)
message("Number of gene combinations: ", ncol(gene_combn))
# Number of gene combinations: 44850

## Find the combination of genes that are co-expressed
message(Sys.time(), " - computing pair_1")
pair_1 <- assay(spe_expr, "expr")[gene_combn[1, ], ]
lobstr::obj_size(pair_1)
# 448.59 MB
message(Sys.time(), " - computing pair_2")
pair_2 <- assay(spe_expr, "expr")[gene_combn[2, ], ]
lobstr::obj_size(pair_2)
# 287.34 MB
message(Sys.time(), " - computing co_expr")
co_expr <- pair_1 & pair_2
message(Sys.time(), " - done computing co_expr")
lobstr::obj_size(co_expr)
# 195.15 MB

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
# 287.39 MB

## Remove objects we don't need anymore
rm(pair_1, pair_2, co_expr)

## Calculate mean of co-expression across groups of interest
group_list <- rafalib::splitit(spe_expr$c2l_combn_unique)
co_expr_means <- do.call(cbind, lapply(group_list, function(ii) {
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
# ENSG00000115944_ENSG00000057757 ENSG00000136193_ENSG00000111912
#                        5.812500                        4.843750
# ENSG00000186462_ENSG00000168734 ENSG00000113013_ENSG00000149091
#                        4.666667                        4.583333
# ENSG00000254772_ENSG00000134243 ENSG00000148798_ENSG00000134982
#                        4.500000                        4.500000

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
        # scale_fill_manual(values = discrete_cell_palette) +
        scale_fill_discrete_qualitative(palette = discrete_cell_palette) +
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
# cols <- unique(co_expr_se$path_groups_colors)
# names(cols) <- sapply(cols, function(x) {
#     names(co_expr_se$path_groups_colors[co_expr_se$path_groups_colors == x])[1]
# })

discrete_cell_palette <- "Dark2"


pdf(
    file.path(dir_plots, paste0("gene_pairs_deconvo_co-expr_", path_name, ".pdf")),
    width = 35, height = 35
)
for(i in seq(0, 90, by = 10)) {
    p <- my_plot_expression(
        sce = co_expr_se,
        genes = names(
            sort(rowData(co_expr_se)$ratio, decreasing = TRUE)
        )[seq_len(10) + i],
        assay = "co_expr",
        ct = "c2l_combn_unique",
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

## Export for cytospace
top_pairs <- names(sort(rowData(co_expr_se)$ratio, decreasing = TRUE))[seq_len(100)]
cyto <- data.frame(
    source_ensembl = gsub("_.*", "", top_pairs),
    source_gene_name = rowData(co_expr_se[top_pairs, ])$gene_name_1,
    target = gsub(".*_", "", top_pairs),
    target_gene_name = rowData(co_expr_se[top_pairs, ])$gene_name_2,
    ratio = rowData(co_expr_se[top_pairs, ])$ratio
)
write.table(cyto, file = file.path(dir_rdata, paste0("cytoscape_co_expr_se_deconvo", path_name, ".txt")), row.names = FALSE, sep = "\t", quote = FALSE)

## Save object for later
saveRDS(
    co_expr_se,
    file = file.path(dir_rdata, paste0("co_expr_se_deconvo", path_name, ".rds"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
