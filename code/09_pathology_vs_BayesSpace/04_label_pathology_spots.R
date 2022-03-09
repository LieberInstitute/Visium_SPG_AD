# sgejobs::job_single(
#     "label_pathology_spots",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "10G",
#     command = "Rscript 04_label_pathology_spots.R",
#     create_logdir = TRUE
# )

library("here")
library("spatialLIBD")
library("sessioninfo")
library("BayesSpace")
library("paletteer")
library("colorblindr") # remotes::install_github("clauswilke/colorblindr")

## Create output directory
dir_rdata <-
    here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
    )
dir_plots <- here::here("plots", "09_pathology_vs_BayesSpace")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE)

## Load the data
spe <- readRDS(
    here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "wholegenome",
        "spe_harmony_wholegenome.rds"
    )
)


## Set pathology levels
spe$path_pTau <-
    ifelse(spe$NpTau > 8 | spe$PpTau > 0.0143, "pTau+", "pTau-")
spe$path_Abeta <-
    ifelse(spe$NAbeta > 1 | spe$PAbeta > 0.108, "Abeta+", "Abeta-")
spe$path_groups <- gsub("au|eta", "", paste0(spe$path_pTau, "_", spe$path_Abeta))
has_path <- which(spe$path_pTau == "pTau+" | spe$path_Abeta == "Abeta+")

## Find neighboring spots
neighbors_list <-
    BayesSpace:::.find_neighbors(spe, platform = "Visium")
## Undo the C++ index from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialCluster.R#L231
path_neighbors_all <- unlist(neighbors_list[has_path]) + 1


which_neighbors <- function(var, values) {
    i <- which(colData(spe)[[var]] %in% values)
    ## Undo the C++ index from
    ## https://github.com/edward130603/BayesSpace/blob/master/R/spatialCluster.R#L231
    res <- sort(unique(unlist(neighbors_list[i]))) + 1

    ## Prioritize the initial pathology categories
    res[!res %in% has_path]
}

path_neighbor_pTau <- which_neighbors("path_pTau", "pTau+")
path_neighbor_Abeta <- which_neighbors("path_Abeta", "Abeta+")
path_neighbor_any <- which_neighbors("path_groups", c("pT+_Ab-", "pT-_Ab+", "pT+_Ab+"))
path_neighbor_both <- intersect(path_neighbor_Abeta, path_neighbor_pTau)

## Compare old code vs new
stopifnot(identical(path_neighbor_any, sort(unique(path_neighbors_all[!path_neighbors_all %in% has_path]))))

## Compare individual pathology results vs any
stopifnot(identical(path_neighbor_any, sort(union(path_neighbor_Abeta, path_neighbor_pTau))))

## Add the neighbors
spe$path_groups[setdiff(path_neighbor_pTau, path_neighbor_both)] <- "next_pT+"
spe$path_groups[setdiff(path_neighbor_Abeta, path_neighbor_both)] <- "next_Ab+"
spe$path_groups[path_neighbor_both] <- "next_both"

## Check that things match
stopifnot(identical(length(path_neighbor_any), length(grep("next_", spe$path_groups))))

## Simplify some groups even more
spe$path_groups <- gsub("pT-_Ab\\+", "Ab+", gsub("pT\\+_Ab-", "pT+", gsub("pT-_Ab-", "none", gsub("pT\\+_Ab\\+", "both", spe$path_groups))))

spe$path_groups <- factor(spe$path_groups, levels = c("none", "Ab+", "next_Ab+", "pT+", "next_pT+", "both", "next_both"))
addmargins(table(spe$path_groups))
#  none       Ab+  next_Ab+       pT+  next_pT+      both next_both       Sum
# 19244      1279      3066      8960      3638       736      1192     38115
round(addmargins(table(spe$path_groups)) / ncol(spe) * 100, 2)
#  none       Ab+  next_Ab+       pT+  next_pT+      both next_both       Sum
# 50.49      3.36      8.04     23.51      9.54      1.93      3.13    100.00

## Load pathology colors
source(here("code", "colors_pathology.R"), echo = TRUE, max.deparse.length = 500)

## Visualize pathology spots
vis_grid_clus(
    spe = spe,
    clustervar = "path_groups",
    pdf_file = file.path(dir_plots, "pathology_groups.pdf"),
    sort_clust = FALSE,
    colors = colors_pathology ,
    spatial = FALSE,
    point_size = 2
)

pdf(file.path(dir_plots, "pathology_groups_colorblind.pdf"), width = 16, height = 14)
p <- vis_clus(
    spe = spe,
    sampleid = "V10T31036_D1_Br3880",
    clustervar = "path_groups",
    colors = colors_pathology ,
    spatial = FALSE,
    point_size = 2
)
colorblindr::cvd_grid(p)
dev.off()

## Export pathology levels for later
for (i in colnames(colData(spe))[grep("^path_", colnames(colData(spe)))]) {
    cluster_export(spe,
        i,
        cluster_dir = dir_rdata
    )
}


## Remove the "both" category in this alternative and prioritize Abeta+ neighbors
## over pTau+
spe$path_groups <- gsub("au|eta", "", paste0(spe$path_pTau, "_", spe$path_Abeta))
table(spe$path_groups)
# pT-_Ab- pT-_Ab+ pT+_Ab- pT+_Ab+
#   27140    1279    8960     736
## Key change number 1: prioritze Abeta+ over pTau+
spe$path_groups[spe$path_groups == "pT+_Ab+"] <- "pT-_Ab+"
table(spe$path_groups)
# pT-_Ab- pT-_Ab+ pT+_Ab-
#   27140    2015    8960

has_path <- which(spe$path_pTau == "pTau+" | spe$path_Abeta == "Abeta+")
## Make has_path into an argument in this version
which_neighbors <- function(var, values, has_path) {
    i <- which(colData(spe)[[var]] %in% values)
    ## Undo the C++ index from
    ## https://github.com/edward130603/BayesSpace/blob/master/R/spatialCluster.R#L231
    res <- sort(unique(unlist(neighbors_list[i]))) + 1

    ## Prioritize the initial pathology categories
    res[!res %in% has_path]
}

path_neighbor_pTau <- which_neighbors("path_pTau", "pTau+", has_path)

## Key change number 2: prioritize Abeta+ neighbors over pTau+
path_neighbor_Abeta <- which_neighbors("path_Abeta", "Abeta+", which(spe$path_Abeta == "Abeta+"))
path_neighbor_both <- intersect(path_neighbor_Abeta, path_neighbor_pTau)

## Add the neighbors
spe$path_groups[setdiff(path_neighbor_pTau, path_neighbor_both)] <- "next_pT+"
spe$path_groups[setdiff(path_neighbor_Abeta, path_neighbor_both)] <- "next_Ab+"
spe$path_groups[path_neighbor_both] <- "next_both"

## Simplify some groups even more
spe$path_groups <- gsub("pT-_Ab\\+", "Ab+", gsub("pT\\+_Ab-", "pT+", gsub("pT-_Ab-", "none", spe$path_groups)))

spe$path_groups <- factor(spe$path_groups, levels = c("none", "Ab+", "next_Ab+", "pT+", "next_pT+", "next_both"))
addmargins(table(spe$path_groups))
#  none       Ab+  next_Ab+       pT+  next_pT+ next_both       Sum
# 19244      2015      5645      6381      3638      1192     38115
round(addmargins(table(spe$path_groups)) / ncol(spe) * 100, 2)
# none       Ab+  next_Ab+       pT+  next_pT+ next_both       Sum
# 50.49      5.29     14.81     16.74      9.54      3.13    100.00


vis_grid_clus(
    spe = spe,
    clustervar = "path_groups",
    pdf_file = file.path(dir_plots, "pathology_groups_alternative.pdf"),
    sort_clust = FALSE,
    colors = colors_pathology[names(colors_pathology) %in% levels(spe$path_groups)],
    spatial = FALSE,
    point_size = 2
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
