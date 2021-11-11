library("spatialLIBD")
library("here")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
# spe <- readRDS("spe_workflow_Visium_spatialLIBD.rds")

# system("ln -s ../../processed-data/spe/spe_raw.Rdata spe_raw.Rdata")
local_path <- here::here("processed-data", "spe", "spe_raw.Rdata")
if(file.exists(local_path)){
    load(local_path, verbose = TRUE)
} else {
    load("spe_raw.Rdata", verbose = TRUE)
}
spe_raw <- spe_raw[, spatialData(spe_raw)$in_tissue]
# rm(spe_raw)

vars <- colnames(colData(spe_raw))

## Deploy the website
spatialLIBD::run_app(
    spe_raw,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "Visium IF AD, Kwon SH et al, 2021",
    spe_discrete_vars = c(
        vars[grep("10x_", vars)],
        # "subject",
        # "sex",
        # "race",
        # "diagnosis",
        # "BCrating",
        # "braak",
        # "cerad",
        # "overlaps_tissue",
        "ManualAnnotation"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        # "age",
        # "pmi",
        # "rin",
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
        "PpTau"
    ),
    default_cluster = "10x_graphclust"
)
