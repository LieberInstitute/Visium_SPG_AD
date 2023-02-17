library("here")
library("sessioninfo")
library("spatialLIBD")

## Locate data directory
dir_rdata <- here::here(
        "code",
        "05_deploy_app_wholegenome"
    )

load(file.path(dir_rdata, "Visium_IF_AD_modeling_results.Rdata"),
    verbose = TRUE
)
sce_pseudo <-
    readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))


## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
