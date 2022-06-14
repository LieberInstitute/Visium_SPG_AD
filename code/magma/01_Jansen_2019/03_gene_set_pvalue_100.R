
## load library
library("spatialLIBD")
library("scater") ## to compute some reduced dimensions
library("dplyr")
library("tidyr")

##Load wholegenome data

dir_rdata <- here::here(
        "processed-data",
        "11_grey_matter_only",
        "without_Br3873")

load(file.path(dir_rdata,
               "wholegenome",
               "Visium_IF_AD_modeling_results.Rdata"), verbose = TRUE)

## extract gene name, pvalues
pvalues <- as_tibble(modeling_results$enrichment) |>
    select(starts_with('p_val') |
<<<<<<< HEAD
    contains('gene'))
=======
    contains('ensembl'))
>>>>>>> 688abd933e142b03c8e937c6ccce53c75088a983

##fix column names
colnames(pvalues) <- gsub(
    "p_value_",
    "",
    colnames(pvalues)
)

colnames(pvalues) <- gsub(
    "+",
    "",
    colnames(pvalues),
    fixed = T)

## apply pivot longer
pvalues <- pvalues|>
<<<<<<< HEAD
    pivot_longer(!gene, names_to = "pathology_type",
=======
    pivot_longer(!ensembl, names_to = "pathology_type",
>>>>>>> 688abd933e142b03c8e937c6ccce53c75088a983
                 values_to = "pvalue")

## group by pathology and filter top 100 smallest values
pvalue_top_100 <- pvalues |>
    arrange(pvalue) |>
    group_by(pathology_type) |>
    slice(1:100)

##rearrange and drop p-value column
pvalue_top_100 <- pvalue_top_100[, c(2,1)]

##rename columns as Set and Gene
colnames(pvalue_top_100) <- c("Set", "Gene")

## write out
<<<<<<< HEAD
write.table(pvalue_top_100, file=here::here("code","12_magma",
=======
write.table(pvalue_top_100, file=here::here("code","magma",
>>>>>>> 688abd933e142b03c8e937c6ccce53c75088a983
                                      "pvalues_top_100.txt"), sep="\t",
            row.names=F, col.names=T, quote=F)
