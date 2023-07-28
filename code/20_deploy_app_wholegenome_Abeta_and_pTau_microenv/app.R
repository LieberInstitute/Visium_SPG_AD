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
    z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.2)))
    write.csv(z[, !grepl("^in_rows", colnames(z))],
        file = file.path(
            dir_rdata,
            "Visium_SPG_AD_wholegenome_model_results_FDR20perc.csv"
        )
    )
    dim(z)
    # [1] 8  13
    table(z$model_type, z$test)
    #                Ab_env Ab_env-both_env Ab_env-none Ab_env-pTau_env both_env-Ab_env none-Ab_env noWM pTau_env-Ab_env
    # anova           0               0           0               0               0           0    1               0
    # enrichment      1               0           0               0               0           0    0               0
    # pairwise        0               1           1               1               1           1    0               1
    z[, !colnames(z) %in% c("in_rows", "results", "in_rows_top20")]
    #        top model_type            test    gene      stat         pval         fdr gene_index     logFC         ensembl
    # 1 9978 enrichment          Ab_env TMEM163 -6.503638 8.038351e-07 0.008020666       1360 -2.869182 ENSG00000152128
    # 2    1   pairwise     none-Ab_env    ISLR  5.341039 1.980015e-05 0.197565867       7045  1.173030 ENSG00000129009
    # 3 9978   pairwise Ab_env-pTau_env TMEM163 -5.395948 1.729706e-05 0.172590047       1360 -2.995160 ENSG00000152128
    # 4 9978   pairwise Ab_env-both_env TMEM163 -5.526484 1.255591e-05 0.125282837       1360 -3.237372 ENSG00000152128
    # 5 9978   pairwise     Ab_env-none    ISLR -5.341039 1.980015e-05 0.197565867       7045  1.173030 ENSG00000129009
    # 6    1   pairwise pTau_env-Ab_env TMEM163  5.395948 1.729706e-05 0.172590047       1360 -2.995160 ENSG00000152128
    # 7    1   pairwise both_env-Ab_env TMEM163  5.526484 1.255591e-05 0.125282837       1360 -3.237372 ENSG00000152128
    # 8    1      anova            noWM TMEM163 14.475536 1.611407e-05 0.160786190       1360        NA ENSG00000152128
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
