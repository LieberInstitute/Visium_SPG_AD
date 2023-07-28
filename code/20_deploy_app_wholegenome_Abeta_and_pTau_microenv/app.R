library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
load("spe.Rdata", verbose = TRUE)
local <- FALSE ## For doing things locally and ability to change input dirs
if (local) {
    ## Local tests
    dir_rdata <- here::here(
        "code",
        "20_deploy_app_wholegenome_Abeta_and_pTau_microenv"
    )
} else {
    ## For shinyapps.io
    dir_rdata <- getwd()
}
load(file.path(dir_rdata, "Visium_SPG_AD_modeling_results.Rdata"),
    verbose = TRUE
)
sce_pseudo <-
    readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))

## Change Braak info based on latest information from LIBD pathology
sce_pseudo$BCrating <- NULL ## This variable was removed from the phenotype table
sce_pseudo$braak <- c("Br3854" = "Stage VI", "Br3873" = "Stage V", "Br3880" = "Stage VI", "Br3874" = "Stage IV")[sce_pseudo$subject]
sce_pseudo$cerad <- c("Br3854" = "Frequent", "Br3873" = "Frequent", "Br3880" = "Frequent", "Br3874" = "None")[sce_pseudo$subject]

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
if (local) {
    fix_csv <- function(df) {
        for (i in seq_len(ncol(df))) {
            if (any(grepl(",", df[, i]))) {
                message(paste(Sys.time(), "fixing column", colnames(df)[i]))
                df[, i] <- gsub(",", ";", df[, i])
            }
        }
        return(df)
    }
    z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
    write.csv(z[, !grepl("^in_rows", colnames(z))],
        file = file.path(
            dir_rdata,
            "Visium_SPG_AD_wholegenome_model_results_FDR5perc.csv"
        )
    )
    dim(z)
    # [1] 318  13
    table(z$model_type, z$test)
    #               Ab_env Ab_env-both Ab_env-none both both-Ab_env none none-Ab_env none-pTau noWM pTau-none
    # anova           0           0           0    0           0    0           0         0    1         0
    # enrichment    275           0           0    9           0    3           0         0    0         0
    # pairwise        0          13           1    0          13    0           1         1    0         1
}

vars <- colnames(colData(spe))
path_vars <- vars[grep("^path_", vars)]
path_vars <- path_vars[!grepl("_colors$", path_vars)]

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "Visium SPG AD (Abeta & pTau microenv), Kwon SH et al, 2023",
    spe_discrete_vars = c(
        path_vars,
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
        "NpTau",
        "PpTau",
        "edge_distance"
    ),
    default_cluster = "path_groups",
    docs_path = "www"
)
