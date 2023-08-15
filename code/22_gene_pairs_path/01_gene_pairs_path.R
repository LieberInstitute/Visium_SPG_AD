library("spatialLIBD")
library("sessioninfo")
spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
sce_pseudo <- spatialLIBD::fetch_data("Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe")

## Subset to genes observed
spe_expr <- spe[rownames(sce_pseudo), spe$diagnosis == "AD"]
assay(spe_expr, "expr") <- counts(spe_expr) > 0








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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
