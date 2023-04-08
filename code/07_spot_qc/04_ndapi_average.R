library(sgejobs)
# sgejobs::job_single(
#     "ndapi_average",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 04_ndapi_average.R",
#     create_logdir = TRUE
# )


## looking for NDAPI outliers, related to https://github.com/LieberInstitute/Visium_SPG_AD/issues/78

## Load remaining required packages
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("sessioninfo")
library("dplyr")

seg_df <- data.frame(
    Percent_Abeta = spe_wholegenome$PAbeta,
    Percent_DAPI = spe_wholegenome$PDAPI,
    Percent_pTau = spe_wholegenome$PpTau,
    Number_Abeta = spe_wholegenome$NAbeta,
    Number_DAPI = spe_wholegenome$NDAPI,
    Number_pTau = spe_wholegenome$NpTau,
    sample_id = spe_wholegenome$sample_id_short,
    check.names = FALSE
)


spe <- readRDS(
    here::here(
        "processed-data", "07_spot_qc",
        paste0("spe_", opt$spetype, "_postqc.rds")
    )
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
