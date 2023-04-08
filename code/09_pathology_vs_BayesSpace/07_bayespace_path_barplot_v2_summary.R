## import required libraries
library("here")
library("sessioninfo")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggrepel")

## output directories
dir_plots <-
    here::here(
        "plots",
        "09_pathology_vs_BayesSpace"
    )
dir.create(dir_plots, showWarnings = FALSE)

suffix <- "wholegenome"

## Load pathology colors
source(
    here("code", "colors_pathology.R"),
    echo = TRUE,
    max.deparse.length = 500
)

spe <- readRDS(here::here(
    "processed-data",
    "07_spot_qc",
    paste0("spe_", suffix, "_postqc.rds")
))

## Import BayesSpace clusters
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

## Import pathology levels
spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
    ),
    prefix = ""
)


## Convert from character to a factor, so they appear in the order
## we want
spe$path_groups <-
    factor(
        spe$path_groups,
        levels = c(
            "none",
            "Ab",
            "n_Ab",
            "pTau",
            "n_pTau",
            "both",
            "n_both"
        )
    )

## Shorten names
spe$sample_id_shorter <- gsub("Br", "", spe$sample_id_short)

## Drop the 3 controls since they are not interesting for this
cluster_df <-
    as.data.frame(colData(spe[, !grepl("3874$", spe$sample_id)]))
bayes_cols <-
    cluster_df |> select(matches("BayesSpace_harmony"))


GM_WM_separate <- cluster_df |>
    select(path_groups, BayesSpace_harmony_k02) |>
    group_by(BayesSpace_harmony_k02) |>
    count(path_groups)

# A tibble: 14 × 3
# Groups:   BayesSpace_harmony_k02 [2]
# BayesSpace_harmony_k02 path_groups     n
# <int> <fct>       <int>
#     1                      1 none         3874
# 2                      1 Ab           1072
# 3                      1 n_Ab         2284
# 4                      1 pTau         8810
# 5                      1 n_pTau       3217
# 6                      1 both          727
# 7                      1 n_both       1102
# 8                      2 none         2461
# 9                      2 Ab            196
# 10                      2 n_Ab          718
# 11                      2 pTau          149
# 12                      2 n_pTau        415
# 13                      2 both            9
# 14                      2 n_both         90

# k = 1 is GM, K =2 is WM

# _____________

## GM and WM grouped together
GMandWM <- cluster_df |>
    select(path_groups, BayesSpace_harmony_k02) |>
    count(path_groups)
# > GMandWM
# path_groups    n
# 1        none 6335
# 2          Ab 1268
# 3        n_Ab 3002
# 4        pTau 8959
# 5      n_pTau 3632
# 6        both  736
# 7      n_both 1192
