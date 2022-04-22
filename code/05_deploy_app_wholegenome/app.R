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
sce_pseudo$APOe <-
    c(
        "Br3854" = "E3/E4",
        "Br3873" = "E3/E3",
        "Br3880" = "E3/E3",
        "Br3874" = "E2/E3"
    )[sce_pseudo$subject]

rowData(sce_pseudo)$gene_search <-
    paste0(rowData(sce_pseudo)$gene_name,
        "; ",
        rowData(sce_pseudo)$gene_id)
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
spe$spatialLIBD <- spe$path_groups
sce_pseudo$layer_guess_reordered_short <- sce_pseudo$path_groups
# pca <- prcomp(t(assays(sce_pseudo)$logcounts))
# reducedDim(sce_pseudo, "PCA") <- pca$x[, seq_len(20)]

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
