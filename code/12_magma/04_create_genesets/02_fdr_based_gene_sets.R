library(sgejobs)
sgejobs::job_single(
    "02_fdr_based_gene_sets",
    create_shell = TRUE,
    queue = "bluejay",
    memory = "20G",
    command = "Rscript 02_fdr_based_gene_sets.R",
    create_logdir = TRUE
)


## load library
library("spatialLIBD")
library("scater") ## to compute some reduced dimensions
library("dplyr")
library("tidyr")
library("sessioninfo")

## Load wholegenome data

dir_rdata <- here::here(
    "processed-data",
    "11_grey_matter_only"
)

load(file.path(
    dir_rdata,
    "wholegenome",
    "Visium_IF_AD_modeling_results.Rdata"
), verbose = TRUE)

## extract gene name, fdr
fdrs <- as_tibble(modeling_results$enrichment) |>
    select(starts_with("fdr") | contains("ensembl"))


## fix column names
colnames(fdrs) <- gsub(
    "fdr_",
    "",
    colnames(fdrs)
)


#############
## apply pivot longer
fdrs_long <- fdrs |>
    pivot_longer(!ensembl,
                 names_to = "pathology_type",
                 values_to = "fdr"
    )

## filter out any value greater than or equal to 0.1 and group by pathology
fdrs_filtered <- fdrs_long |>
    filter(fdr< 0.1) |>
    group_by(pathology_type)

nrow(fdrs_filtered |> filter(pathology_type == "n_Ab"))
nrow(fdrs_filtered |> filter(pathology_type == "Ab"))

## rearrange and drop p-value column
fdr_gene_set <- fdrs_filtered[, c(2, 1)]

## rename columns as Set and Gene
colnames(fdr_gene_set) <- c("Set", "Gene")

## write out
write.table(fdr_gene_set,
            file = here::here(
                "code", "12_magma",
                "fdr_gene_set.txt"
            ), sep = "\t",
            row.names = F, col.names = T, quote = F
)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
