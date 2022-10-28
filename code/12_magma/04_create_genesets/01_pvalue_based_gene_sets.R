# library(sgejobs)
# sgejobs::job_single(
#     "pvalue_based_gene_sets",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 01_pvalue_based_gene_sets.R",
#     create_logdir = TRUE
# )


## load library
library("spatialLIBD")
library("scater") ## to compute some reduced dimensions
library("dplyr")
library("tidyr")

## Load wholegenome data

dir_rdata <- here::here(
    "processed-data",
    "11_grey_matter_only",
    "without_Br3873"
)

load(file.path(
    dir_rdata,
    "wholegenome",
    "Visium_IF_AD_modeling_results.Rdata"
), verbose = TRUE)

## extract gene name, pvalues
pvalues <- as_tibble(modeling_results$enrichment) |>
    select(starts_with("p_value") | contains("ensembl"))


## fix column names
colnames(pvalues) <- gsub(
    "p_value_",
    "",
    colnames(pvalues)
)

colnames(pvalues) <- gsub(
    "+",
    "",
    colnames(pvalues),
    fixed = T
)
#############
## apply pivot longer
pvalues <- pvalues |>
    pivot_longer(!ensembl,
        names_to = "pathology_type",
        values_to = "pvalue"
    )

## group by pathology and filter top 100 smallest values
pvalue_top_100 <- pvalues |>
    arrange(pvalue) |>
    group_by(pathology_type) |>
    slice(1:100)

## rearrange and drop p-value column
pvalue_top_100 <- pvalue_top_100[, c(2, 1)]

## rename columns as Set and Gene
colnames(pvalue_top_100) <- c("Set", "Gene")

## write out
write.table(pvalue_top_100,
    file = here::here(
        "code", "magma",
        "pvalues_top_100.txt"
    ), sep = "\t",
    row.names = F, col.names = T, quote = F
)

################## 50 genes##################
pvalue_top_50 <- pvalues |>
    arrange(pvalue) |>
    group_by(pathology_type) |>
    slice(1:50)

## rearrange and drop p-value column
pvalue_top_50 <- pvalue_top_50[, c(2, 1)]

## rename columns as Set and Gene
colnames(pvalue_top_50) <- c("Set", "Gene")

## write out
write.table(pvalue_top_50,
    file = here::here(
        "code", "magma",
        "pvalues_top_50.txt"
    ), sep = "\t",
    row.names = F, col.names = T, quote = F
)

################ 200 genes###################
pvalue_top_200 <- pvalues |>
    arrange(pvalue) |>
    group_by(pathology_type) |>
    slice(1:200)

## rearrange and drop p-value column
pvalue_top_200 <- pvalue_top_200[, c(2, 1)]

## rename columns as Set and Gene
colnames(pvalue_top_200) <- c("Set", "Gene")

## write out
write.table(pvalue_top_200,
    file = here::here(
        "code", "magma",
        "pvalues_top_200.txt"
    ), sep = "\t",
    row.names = F, col.names = T, quote = F
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
