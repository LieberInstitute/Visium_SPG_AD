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
load("spe_targeted.Rdata", verbose = TRUE)

vars <- colnames(colData(spe_targeted))

## Deploy the website
spatialLIBD::run_app(
    spe_targeted,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "Visium IF AD (TGE), Kwon SH et al, 2021",
    spe_discrete_vars = c(
        vars[grep("10x_", vars)],
        "ManualAnnotation",
        "edge_spots",
        vars[grep("^scran_", vars)]
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
    default_cluster = "10x_graphclust"
)
