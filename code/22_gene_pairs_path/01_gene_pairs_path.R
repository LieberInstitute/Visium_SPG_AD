library("spatialLIBD")
spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
sce_pseudo <- spatialLIBD::fetch_data("Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe")

## Subset to genes observed
spe_expr <- spe[rownames(sce_pseudo), spe$diagnosis == "AD"]
assay(spe_expr, "expr") <- counts(spe_expr) > 0

gene_combn <- combn(10, 2)
ncol(gene_combn)

Sys.time()
co_expr_loop <- do.call(rbind, apply(gene_combn, 2, function(x) {
    pair <- assay(spe_expr, "expr")[x, ]
    pair[1, , drop = FALSE] & pair[2, , drop = FALSE]
}))
Sys.time()


pair_1 <- assay(spe_expr, "expr")[gene_combn[1, ], ]
pair_2 <- assay(spe_expr, "expr")[gene_combn[2, ], ]
co_expr_v2 <- pair_1 & pair_2
identical(co_expr_loop, co_expr_v2)






gene_combn <- combn(500, 2)
ncol(gene_combn)
# [1] 124750

Sys.time()
# [1] "2023-08-15 12:32:14 EDT"
pair_1 <- assay(spe_expr, "expr")[gene_combn[1, ], ]
Sys.time()
# [1] "2023-08-15 12:32:15 EDT"
lobstr::obj_size(pair_1)
# 2.91 GB
Sys.time()
# [1] "2023-08-15 12:32:15 EDT"
pair_2 <- assay(spe_expr, "expr")[gene_combn[2, ], ]
Sys.time()
# [1] "2023-08-15 12:32:22 EDT"
lobstr::obj_size(pair_2)
# 2.71 GB
Sys.time()
# [1] "2023-08-15 12:32:22 EDT"
co_expr <- pair_1 & pair_2
# Error: vector memory exhausted (limit reached?)
Sys.time()
# [1] "2023-08-15 12:32:56 EDT"
lobstr::obj_size(co_expr)
# 389.46 MB
