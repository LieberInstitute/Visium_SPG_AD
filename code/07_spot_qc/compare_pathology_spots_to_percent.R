## import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")

## Load basic SPE data
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)
load(here::here("processed-data", "07_spot_qc", "spe_targeted_postqc.Rdata"), verbose = TRUE)



## output directories
dir_plots <- here::here("plots", "07_spot_qc")


## Create function to create 2x2 Table
two_by_two_table <- function(spe_object, sample_num, pathology, n = 0, p = 0.01) {
    # n = threshold for number of pathology 'blobs' in spot
    # p = percentage of pathology pixels in spot
    # add optional params (diagnosis)

    path_df <- data.frame(
        spot_id = rownames(colData(spe_object)),
        diagnosis = colData(spe_object)$diagnosis,
        sample_id = colData(spe_object)$sample_id,
        NAbeta = colData(spe_object)$NAbeta,
        NpTau = colData(spe_object)$NpTau,
        PAbeta = colData(spe_object)$PAbeta,
        PpTau = colData(spe_object)$PpTau
    )

    path_df <- path_df |> filter(sample_id == sample_num)

    ## 2x2 table column 1: Number present, column 2: Number absent
    ## row 1: Percent > 0.01, column 2: Percent absent

    if (pathology == "Abeta") {
        col1_row1 <- path_df |>
            filter(NAbeta > n & PAbeta > p) |>
            count()
        col1_row2 <- path_df |>
            filter(NAbeta > n & PAbeta <= p) |>
            count()
        col2_row1 <- path_df |>
            filter(NAbeta == n & PAbeta > p) |>
            count()
        col2_row2 <- path_df |>
            filter(NAbeta == n & PAbeta <= p) |>
            count()

        row_names <- c("PAbeta Present", "PAbeta Absent")
        col_names <- c("NAbeta Present", "NAbeta Absent")
    }

    if (pathology == "pTau") {
        col1_row1 <- path_df |>
            filter(NpTau > n & PpTau > p) |>
            count()
        col1_row2 <- path_df |>
            filter(NpTau > n & PpTau <= p) |>
            count()
        col2_row1 <- path_df |>
            filter(NpTau == n & PpTau > p) |>
            count()
        col2_row2 <- path_df |>
            filter(NpTau == n & PpTau <= p) |>
            count()

        row_names <- c("PpTau Present", "PpTau Absent")
        col_names <- c("NpTau Present", "NpTau Absent")
    }

    vals <- c(col1_row1, col1_row2, col2_row1, col2_row2)
    matrix_2x2 <- matrix(vals, nrow = 2, ncol = 2, dimnames = list(row_names, col_names))
    return(matrix_2x2)
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
