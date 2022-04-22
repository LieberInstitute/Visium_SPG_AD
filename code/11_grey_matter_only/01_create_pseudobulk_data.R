# library(sgejobs)

# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "01_create_pseudobulk_data",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 01_create_pseudobulk_data.sh

# Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## For testing
if (FALSE) {
    opt <- list(spetype = "wholegenome")
}


library("SpatialExperiment")
library("here")
library("spatialLIBD")
library("scuttle")
library("edgeR")
library("sessioninfo")
library("jaffelab")

## output directory
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

## load spe data
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            opt$spetype,
            paste0("spe_harmony_", opt$spetype, ".rds")
        )
    )

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        opt$spetype,
        "clusters_BayesSpace"
    ),
    prefix = ""
)

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
    ),
    prefix = ""
)

## subset spe data based on subject and cluster 1 for k = 2
spe_new <- spe[, !spe$subject == "Br3874"]

if (opt$spetype == "wholegenome") {
    spe_new <- spe_new[, spe_new$BayesSpace_harmony_k02 != 2]
} else {
    spe_new <- spe_new[, spe_new$BayesSpace_harmony_k04 != 4]
}

# > dim(colData(spe))
# [1] 38115   111
# # > dim(colData(spe_new))
# [1] 21086   111

## > unique(spe$path_groups)
# [1] "none"      "next_Ab+"  "Ab+"       "both"      "pT+"       "next_both"
# [7] "next_pT+"


## pseudobulk across pathology labels
sce_pseudo <- aggregateAcrossCells(
    spe_new,
    DataFrame(
        path_groups = spe_new$path_groups,
        sample_id = spe_new$sample_id
    )
)
x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)

stopifnot(identical(rownames(x), rownames(sce_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(sce_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <- x

## From
## https://github.com/LieberInstitute/spatialDLPFC/blob/e38213e47f780074af6a4575b404765a486590e6/code/analysis/09_region_differential_expression/preliminary_analysis.R#L47-L55
rowData(sce_pseudo)$low_expr <- filterByExpr(sce_pseudo)
summary(rowData(sce_pseudo)$low_expr)
sce_pseudo <- sce_pseudo[which(!rowData(sce_pseudo)$low_expr), ]
dim(sce_pseudo)

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68
pca <- prcomp(t(assays(sce_pseudo)$logcounts))
message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(sce_pseudo)
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

## We don't want to model the pathology groups as integers / numeric
## so let's double check this
stopifnot(is.factor(sce_pseudo$path_groups) || is.character(sce_pseudo$path_groups))

## Add APOe genotype info
sce_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E3", "Br3880" = "E3/E3", "Br3874" = "E2/E3")[sce_pseudo$subject]

## For the spatialLIBD shiny app
rowData(sce_pseudo)$gene_search <-
    paste0(
        rowData(sce_pseudo)$gene_name,
        "; ",
        rowData(sce_pseudo)$gene_id
    )
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
## Something I need to fix on the shiny app since it's hardcoded to use
## this variable right now
sce_pseudo$layer_guess_reordered_short <- sce_pseudo$path_groups

## save RDS file
saveRDS(
    sce_pseudo,
    file = here::here(
        "processed-data", "11_grey_matter_only", opt$spetype,
        paste0("sce_pseudo_pathology_", opt$spetype, ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
