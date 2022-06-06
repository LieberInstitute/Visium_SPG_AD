library("spatialLIBD")
library("scater") ## to compute some reduced dimensions
library("dplyr")

#### Load the data ####

dir_rdata <- here::here(
        "processed-data",
        "11_grey_matter_only",
        "without_Br3873",
        "wholegenome")

load(file.path(dir_rdata, "Visium_IF_AD_modeling_results.Rdata"), verbose = TRUE)
sce_pseudo <- readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))

#### Fix column names ####
colnames(modeling_results$enrichment) <- gsub(
    "pos",
    "+",
    colnames(modeling_results$enrichment)
)

## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

enrichment_results <- modeling_results$enrichment


