library("spatialLIBD")
library("scater") ## to compute some reduced dimensions
library("dplyr")
library("tidyr")

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


#### extract gene name, pvalues ####
pvalues <- as_tibble(modeling_results$enrichment) |>
    select(starts_with('p_val') |
    contains('gene'))

colnames(pvalues) <- gsub(
    "p_value_",
    "",
    colnames(pvalues)
)

#### apply pivot longer ####
pvalues <- pvalues|>
    pivot_longer(!gene, names_to = "pathology_type",
                 values_to = "pvalue")



