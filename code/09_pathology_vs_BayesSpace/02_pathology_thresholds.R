## import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("tidyr")

## Load basic SPE data
spe <- readRDS(
    here::here(
        "processed-data", "07_spot_qc", "spe_wholegenome_postqc.rds"
    )
)

controls <- c("V10A27106_A1_Br3874", "V10A27004_A1_Br3874", "V10T31036_A1_Br3874")
# last one shouldn't be used for pTau


## find max of NpTau, PpTau

path_df <- data.frame(
    spot_id = colnames(spe),
    diagnosis = spe$diagnosis,
    sample_id = spe$sample_id,
    NAbeta = spe$NAbeta,
    NpTau = spe$NpTau,
    PAbeta = spe$PAbeta,
    PpTau = spe$PpTau
)

## Just for NpTau/PpTau
path_df |>
    dplyr::filter(sample_id %in% controls) |>
    summarise_if(is.numeric, max, na.rm = TRUE)

#   NAbeta NpTau    PAbeta      PpTau
# 1      4     8 0.1983471 0.01433482

## Just for NAbeta/PAbeta
path_df |>
    dplyr::filter(sample_id %in% controls[c(1, 3)]) |>
    summarise_if(is.numeric, max, na.rm = TRUE)

#   NAbeta NpTau   PAbeta      PpTau
# 1      3     8 0.149126 0.01433482


## Frequency of unique NAbeta values across all controls

path_df |>
    dplyr::filter(sample_id %in% controls) |>
    count(NAbeta) |>
    group_by(NAbeta) |>
    mutate(prop = prop.table(n))

# '''
#     NAbeta   n  prop
#    <int> <int> <dbl>
# 1      0 12963     1
# 2      1    22     1
# 3      2     3     1
# 4      3     2     1
# 5      4     1     1
# '''


## Quantiles for NAbeta
path_df |>
    dplyr::filter(sample_id %in% controls) |>
    group_by(sample_id) |>
    summarise(q = list(quantile(NAbeta)), na.rm = TRUE) |>
    unnest_wider(q)

# '''
# sample_id            `0%` `25%` `50%` `75%` `100%` na.rm
#   <chr>               <dbl> <dbl> <dbl> <dbl>  <dbl> <lgl>
# 1 V10A27004_A1_Br3874     0     0     0     0      4 TRUE
# 2 V10A27106_A1_Br3874     0     0     0     0      2 TRUE
# 3 V10T31036_A1_Br3874     0     0     0     0      3 TRUE
# '''


## New percentiles for NAbeta
path_df |>
    dplyr::filter(sample_id %in% controls) |>
    group_by(sample_id) |>
    summarise(
        percentiles = scales::percent(c(0.95, 0.96, 0.97, 0.98, 0.99, 0.999)),
        NAbeta = quantile(NAbeta, c(0.95, 0.96, 0.97, 0.98, 0.99, 0.999)),
        na.rm = TRUE
    )
# Everything zero except 0.999 where NAbeta = 1






## Quantiles for PAbeta
path_df |>
    dplyr::filter(sample_id %in% controls) |>
    group_by(sample_id) |>
    summarise(q = list(quantile(PAbeta)), na.rm = TRUE) |>
    unnest_wider(q)

# '''
# sample_id            `0%` `25%` `50%` `75%` `100%` na.rm
# <chr>               <dbl> <dbl> <dbl> <dbl>  <dbl> <lgl>
# 1 V10A27004_A1_Br3874     0     0     0     0 0.198  TRUE
# 2 V10A27106_A1_Br3874     0     0     0     0 0.0649 TRUE
# 3 V10T31036_A1_Br3874     0     0     0     0 0.149  TRUE
# '''

## New percentiles for PAbeta
path_df |>
    dplyr::filter(sample_id %in% controls) |>
    group_by(sample_id) |>
    summarise(
        percentiles = scales::percent(c(0.95, 0.96, 0.97, 0.98, 0.99, 0.999)),
        NAbeta = quantile(PAbeta, c(0.95, 0.96, 0.97, 0.98, 0.99, 0.999)),
        na.rm = TRUE
    )
## for 004 and 1036 99.9% is 0.108 and 0.0543 respectively. Zeros for everything else.


path_df_AD <- path_df |> dplyr::filter(!sample_id %in% controls)
count(path_df_AD) # 25124 total spots in all AD samples

thresholded <- path_df_AD |> dplyr::filter(NAbeta > 1 | PAbeta > 0.108)
count(thresholded)
# 1 2004

path_df_AD |>
    dplyr::filter(NAbeta >= 1 | PAbeta >= 0.108) |>
    count()
#      n
# 1 2861

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
