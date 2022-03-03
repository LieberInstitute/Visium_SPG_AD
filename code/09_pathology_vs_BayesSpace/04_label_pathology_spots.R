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

## Create output directory
dir_rdata <- here::here("processed-data", "09_pathology_vs_BayesSpace", "pathology_levels")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load the data
spe <- readRDS(
    here::here(
        "processed-data", "07_spot_qc", "spe_wholegenome_postqc.rds"
    )
)

## Set pathology levels
spe$path_pTau <- ifelse(spe$NpTau > 7 | spe$PpTau > 0.014, "pTau+", "pTau-")
spe$path_Abeta <- ifelse(spe$NAbeta > 1 | spe$PAbeta > 0.108, "Abeta+", "Abeta-")
spe$path_groups <- paste0(spe$path_pTau, "_", spe$path_Abeta)

## TODO find neighbors

## Export pathology levels for later
for (i in colnames(colData(spe))[grep("^path_", colnames(colData(spe)))]) {
    cluster_export(
        spe,
        i,
        cluster_dir = dir_rdata
    )
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
