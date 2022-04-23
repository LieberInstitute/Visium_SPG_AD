library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
load("spe.Rdata", verbose = TRUE)
load("Visium_IF_AD_modeling_results.Rdata", verbose = TRUE)
sce_pseudo <- readRDS("sce_pseudo_pathology_wholegenome.rds")

## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

## Extract FDR < 5%
## From
## https://github.com/LieberInstitute/brainseq_phase2/blob/be2b7f972bb2a0ede320633bf06abe1d4ef2c067/supp_tabs/create_supp_tables.R#L173-L181
# fix_csv <- function(df) {
#     for (i in seq_len(ncol(df))) {
#         if (any(grepl(",", df[, i]))) {
#             message(paste(Sys.time(), "fixing column", colnames(df)[i]))
#             df[, i] <- gsub(",", ";", df[, i])
#         }
#     }
#     return(df)
# }
# z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
# write.csv(z, file = "Visium_IF_AD_wholegenome_model_results_FDR5perc.csv")

vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "Visium IF AD, Kwon SH et al, 2022",
    spe_discrete_vars = c(
        vars[grep("^path_", vars)],
        "ManualAnnotation",
        vars[grep("^BayesSpace_", vars)],
        vars[grep("^graph_", vars)],
        "edge_spots",
        vars[grep("^scran_", vars)],
        vars[grep("^10x_", vars)]
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "NAbeta",
        "PAbeta",
        "NDAPI",
        "PDAPI",
        "NGFAP",
        "PGFAP",
        "NLipofuscin",
        "PLipofuscin",
        "NMAP2",
        "PMAP2",
        "NpTau",
        "PpTau",
        "edge_distance"
    ),
    default_cluster = "path_groups"
)
