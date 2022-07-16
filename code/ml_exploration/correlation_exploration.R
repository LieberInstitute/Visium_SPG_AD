# library(sgejobs)
# sgejobs::job_single(
#         "correlation_exploration",
#         create_shell = TRUE,
#         queue = "bluejay",
#         memory = "40G",
#         command = "Rscript correlation_exploration.R",
#         create_logdir = TRUE
#     )

library(rtracklayer)
library(GenomicRanges)
library(SpatialExperiment)
library(here)
library(dplyr)
library(Matrix)


# library(usethis)
# usethis::edit_r_environ()


## first work with correlation plots
spe_wholegenome <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            "wholegenome",
            paste0("spe_harmony_", "wholegenome", ".rds")
        )
    )


# rm(spe_counts)
#
#
# spe_counts <- counts(spe_wholegenome)
# spe_counts <- as.matrix(spe_counts)
#
# spe_colData <- t(as.matrix(colData(spe_wholegenome)))
#
# spe_colData  <- spe_colData[rownames(spe_colData) %in%
#                          c("sample_id_short", "NAbeta", "PAbeta",
#                            "NpTau", "PpTau"), ]
# rm(spe_wholegenome)
# spe_counts <- rbind(spe_counts, spe_colData)
#
#
# corrs <- vector(mode = "list", length = 27853)
#
# for (i in nrow(spe_counts)){
#     append(corrs,
#            cor( spe_counts[i, ,drop =FALSE],
#                 as.numeric(spe_colData["NAbeta"])))
# }
#
# cor( spe_counts[1,1:38115,drop =FALSE],
#      spe_colData$NAbeta)
# # > dim(spe_colData)
# # [1]     5 38115
#
# #
# # rm(spe_wholegenome)
# #
# # rbind(spe_counts,spe_colData)
#
#



# Object obj1 is the Seurat object having the highest number of cells
# Object obj2 is the second Seurat object with lower number of cells

# Compute the length of cells from obj2
# cells.to.sample <- seq(1:38115)
#
# # Sample from obj1 as many cells as there are cells in obj2
# # For reproducibility, set a random seed
# set.seed(111)
# sampled.cells <- sample(x = cells.to.sample ,
#                         size = 1000 ,
#                         replace = F)

spe_counts <- counts(spe_wholegenome)
# > dim(spe_counts)
# [1] 27853  1000

spe_cols <- colData(spe_wholegenome)
# > dim(spe_cols)
# [1] 1000   81

rm(spe_wholegenome)
spe_cols_mat <- t(as.matrix(spe_cols))

spe_counts <- as.matrix(spe_counts)

spe_counts <- rbind(spe_counts, spe_cols_mat)

dir.create(here("processed-data", "ml_exploration"))
saveRDS(spe_counts, file = here("processed-data", "ml_exploration", "counts_and_col.RDS"))
