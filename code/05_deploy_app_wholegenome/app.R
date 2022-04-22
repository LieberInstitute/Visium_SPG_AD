library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
load("spe.Rdata", verbose = TRUE)
load("Visium_IF_AD_modeling_results.Rdata", verbose = TRUE)
sce_pseudo <- readRDS("sce_pseudo_pathology_wholegenome.rds")

spe$spatialLIBD <- spe$path_groups
vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes_extract_all(
        n = nrow(sce_pseudo),
        modeling_results = modeling_results,
        sce_layer = sce_pseudo
    ),
    title = "Visium IF AD, Kwon SH et al, 2021",
    spe_discrete_vars = c(
        vars[grep("^path_", vars)],
        "ManualAnnotation",
        vars[grep("^BayesSpace_", vars)],
        vars[grep("^graph_", vars)],
        "edge_spots",
        vars[grep("^scran_", vars)],
        vars[grep("^10x_", vars)],
        "spatialLIBD"
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
