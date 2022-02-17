##import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("tidyr")
## Load basic SPE data
load(here::here("processed-data", "07_spot_qc", "spe_postqc.Rdata"), verbose = TRUE)




controls <- c("V10A27106_A1_Br3874", "V10A27004_A1_Br3874", "V10T31036_A1_Br3874")
#last one shouldn't be used for pTau


##find max of NpTau, PpTau

path_df <- data.frame(
    spot_id =rownames(colData(spe)),
    diagnosis = colData(spe)$diagnosis,
    sample_id = colData(spe)$sample_id,
    NAbeta = colData(spe)$NAbeta,
    NpTau = colData(spe)$NpTau,
    PAbeta = colData(spe)$PAbeta,
    PpTau = colData(spe)$PpTau)


path_df |> filter(sample_id %in% controls[1:2]) |> summarise_if(is.numeric, max, na.rm = TRUE)

#NAbeta NpTau    PAbeta       PpTau
#  4     7     0.1983471      0.01396914



##Quantiles for NAbeta
path_df |> filter(sample_id %in% controls) |> group_by(sample_id) |>
    summarise( q = list(quantile(NAbeta)), na.rm = TRUE) |> unnest_wider(q)

'''
sample_id            `0%` `25%` `50%` `75%` `100%` na.rm
  <chr>               <dbl> <dbl> <dbl> <dbl>  <dbl> <lgl>
1 V10A27004_A1_Br3874     0     0     0     0      4 TRUE
2 V10A27106_A1_Br3874     0     0     0     0      2 TRUE
3 V10T31036_A1_Br3874     0     0     0     0      3 TRUE
'''

##Frequency of unique NAbeta values across all controls

path_df |> filter(sample_id %in% controls) |> count(NAbeta) |>
    group_by(NAbeta) |> mutate(prop = prop.table(n))

'''
    NAbeta   n  prop
   <int> <int> <dbl>
1      0 12963     1
2      1    22     1
3      2     3     1
4      3     2     1
5      4     1     1
'''

## Quantiles for PAbeta
path_df |> filter(sample_id %in% controls) |> group_by(sample_id) |>
    summarise( q = list(quantile(PAbeta)), na.rm = TRUE) |> unnest_wider(q)

'''
sample_id            `0%` `25%` `50%` `75%` `100%` na.rm
<chr>               <dbl> <dbl> <dbl> <dbl>  <dbl> <lgl>
1 V10A27004_A1_Br3874     0     0     0     0 0.198  TRUE
2 V10A27106_A1_Br3874     0     0     0     0 0.0649 TRUE
3 V10T31036_A1_Br3874     0     0     0     0 0.149  TRUE
'''


