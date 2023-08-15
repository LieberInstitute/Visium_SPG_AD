library("spatialLIBD")
library("lobstr")
library("rafalib")
library("sessioninfo")

## Read in the data
spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
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
# rm(spe)

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

## Remove objects we don't need anymore
rm(pair_1, pair_2)

## Calculate mean of co-expression across groups of interest
path_list <- rafalib::splitit(spe_expr$path_groups)
co_expr_means <- do.call(cbind, lapply(path_list, function(ii) {
    rowMeans(co_expr[, ii])
}))
dim(co_expr_means)

## Find highest (max) and second highest for each gene pair
highest <- apply(co_expr_means, 1, max)
second <- apply(co_expr_means, 1, function(x) {
    max(x[x != max(x)])
})
ratio <- highest / second
message("Summary of ratio of highest / second highest:")
summary(ratio)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
