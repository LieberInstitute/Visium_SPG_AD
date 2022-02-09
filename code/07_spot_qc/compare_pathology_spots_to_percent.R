#import required libraries
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



n_threshold = 0   #number threshold
p_threshold = 0.01  #percentage threshold


##pathology_df


''' path_df_spe <- data.frame(
    spot_id =rownames(colData(spe)),
    sample_id = colData(spe)$sample_id,
    NAbeta = colData(spe)$NAbeta,
    NpTau = colData(spe)$NpTau,
    PAbeta = colData(spe)$PAbeta,
    PpTau = colData(spe)$PpTau
    )
'''



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
