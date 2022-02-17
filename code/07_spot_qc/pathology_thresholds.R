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



