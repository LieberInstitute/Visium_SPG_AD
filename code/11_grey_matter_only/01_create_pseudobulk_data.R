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
check_code <- FALSE
if (FALSE) {
    opt <- list(spetype = "wholegenome")
    check_code <- TRUE
}


library("SpatialExperiment")
library("here")
library("spatialLIBD")
library("scuttle")
library("edgeR")
library("sessioninfo")
library("jaffelab")
library("scater")

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

## Convert from character to a factor, so they appear in the order
## we want
spe$path_groups <-
    factor(
        spe$path_groups,
        levels = c(
            "none",
            "Ab+",
            "next_Ab+",
            "pT+",
            "next_pT+",
            "both",
            "next_both"
        )
    )

## subset spe data based on subject and cluster 1 for k = 2
spe <- spe[, !spe$subject %in% c("Br3874")]

if (opt$spetype == "wholegenome") {
    spe <- spe[, spe$BayesSpace_harmony_k02 != 2]
} else {
    spe <- spe[, spe$BayesSpace_harmony_k04 != 4]
}
dim(spe)

## > unique(spe$path_groups)
# [1] "none"      "next_Ab+"  "Ab+"       "both"      "pT+"       "next_both"
# [7] "next_pT+"
summary(spe$path_groups)
stopifnot(is.factor(spe$path_groups))


## pseudobulk across pathology labels
sce_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        path_groups = spe$path_groups,
        sample_id = spe$sample_id
    )
)
dim(sce_pseudo)

if (check_code) {
    ## Check that we get the same results vs our manual pseudo-bulk code
    ## Original pseudo-bulking code:
    ## https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L89-L162

    library("rafalib")
    pathIndexes <- splitit(paste0(spe$sample_id, '_', spe$path_groups))

    umiComb <-sapply(pathIndexes, function(ii)
            rowSums(assays(spe)$counts[, ii, drop = FALSE]))
    stopifnot(identical(dim(umiComb), dim(sce_pseudo)))

    library("jaffelab")
    path_df <- data.frame(
        sample_id = paste0(ss(colnames(umiComb), '_', 1), "_", ss(colnames(umiComb), '_', 2), "_", ss(colnames(umiComb), '_', 3)),
        path_groups = factor(ss(colnames(umiComb), 'Br[0-9]+_', 2), levels = levels(spe$path_groups)),
        stringsAsFactors = FALSE
    )
    m_layer <- match(path_df$sample_id, spe$sample_id)
    path_df$subject <- spe$subject[m_layer]
    rownames(path_df) <- colnames(umiComb)

    sce_path <- SingleCellExperiment(
        list(counts = umiComb),
        colData = path_df,
        rowData = rowData(spe)
    )

    ## Objects are not in the same order
    stopifnot(identical(colnames(umiComb), colnames(sce_path)))
    m_sce <- match(
        paste0(sce_pseudo$sample_id, "_", sce_pseudo$path_groups),
        colnames(sce_path)
    )
    stopifnot(!any(is.na(m_sce)))
    sce_path <- sce_path[, m_sce]

    ## Now they are in the same order
    stopifnot(identical(sce_pseudo$sample_id, sce_path$sample_id))
    colnames(sce_pseudo) <- colnames(sce_path) ## sce_pseudo has no colnames
    ## making them identical helps with other checks later on

    ## Excellent, the counts are the same!
    stopifnot(identical(counts(sce_pseudo), counts(sce_path)))

    ## Filter by ncells here too
    sce_path <- sce_path[, sce_pseudo$ncells >= 10]
}

## Drop combinations that are very low (very few spots were pseudo-bulked)
## From
## http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#performing-the-de-analysis
summary(sce_pseudo$ncells)
table(sce_pseudo$ncells >= 10)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]

## From
## https://github.com/LieberInstitute/spatialDLPFC/blob/e38213e47f780074af6a4575b404765a486590e6/code/analysis/09_region_differential_expression/preliminary_analysis.R#L47-L55
rowData(sce_pseudo)$low_expr <- filterByExpr(sce_pseudo)
rowData(sce_pseudo)$low_expr_group_sample_id <- filterByExpr(sce_pseudo, group = sce_pseudo$sample_id)
rowData(sce_pseudo)$low_expr_group_path_groups <- filterByExpr(sce_pseudo, group = sce_pseudo$path_groups)
with(rowData(sce_pseudo), table(low_expr, low_expr_group_path_groups))
with(rowData(sce_pseudo), table(low_expr_group_sample_id, low_expr_group_path_groups))

## Check what's the issue with SNAP25
if (check_code) {
    i <- which(rowData(sce_pseudo)$gene_name == "SNAP25")
    summary(counts(sce_pseudo)[i, ])
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # 24.0   231.5   750.0  2001.0  2542.5 10392.0

    rowData(sce_pseudo)[i, ]
    # DataFrame with 1 row and 10 columns
    # source     type         gene_id gene_version   gene_name      gene_type            gene_search  low_expr low_expr_group_sample_id
    # <factor> <factor>     <character>  <character> <character>    <character>            <character> <logical>                <logical>
    #     ENSG00000132639   HAVANA     gene ENSG00000132639           12      SNAP25 protein_coding SNAP25; ENSG00000132..      TRUE                     TRUE
    # low_expr_group_path_groups
    # <logical>
    #     ENSG00000132639                       TRUE

    ?filterByExpr
    # Value
    # Logical vector of length nrow(y) indicating which rows of y to keep in the analysis.

    ## We were originally using !filterByExpr() which meant that we dropped all
    ## highly expressed genes and kept only the ones with low expression values!
    ## see
    ## https://github.com/LieberInstitute/Visium_IF_AD/blob/3c7a829dc2bee97d8d7850c9df40e2eec6335ef0/code/11_grey_matter_only/01_create_pseudobulk_data.R#L142
}

## Check it's on the same order
stopifnot(identical(rownames(sce_pseudo), names(rowData(sce_pseudo)$low_expr_group_path_groups)))
if (check_code) {
    sce_path <- sce_path[rowData(sce_pseudo)$low_expr_group_path_groups, ]
    dim(sce_path)
}

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$low_expr_group_path_groups, ]
dim(sce_pseudo)

## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <-
    edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
        log = TRUE,
        prior.count = 1)

if (check_code) {
    sce_path <- logNormCounts(sce_path)

    ## We used to have large differences in the log counts
    ## related to the !filterByExpr() vs filterByExpr() issue noted earlier
    ## from
    ## https://github.com/LieberInstitute/Visium_IF_AD/blob/3c7a829dc2bee97d8d7850c9df40e2eec6335ef0/code/11_grey_matter_only/01_create_pseudobulk_data.R#L142
    summary(logcounts(sce_path)[, 1] - logcounts(sce_pseudo)[, 1])
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -5.474  -5.474  -5.474  -5.424  -5.338  -5.279
    #
    # ## Updated info with filterByExpr()
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4439 -0.4419 -0.4399 -0.4379 -0.4367 -0.3111

    for (i in seq_len(ncol(sce_path))) {
        print(summary(logcounts(sce_path)[, i] - logcounts(sce_pseudo)[, i]))
    }
    ## There are some differences still
    ## but not as big as they used to be (they used to be at like -5)
    ## So I feel comfortable proceeding with this version
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4439 -0.4419 -0.4399 -0.4379 -0.4367 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3355 -0.3351 -0.3348 -0.3346 -0.3343 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2461 -0.2448 -0.2453 -0.2438 -0.2429
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3065 -0.3064 -0.3065 -0.3063 -0.3063
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3031 -0.3029 -0.3032 -0.3028 -0.3027
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2753 -0.2745 -0.2748 -0.2740 -0.2735
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3051 -0.3050 -0.3051 -0.3049 -0.3048
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3482 -0.3476 -0.3471 -0.3429 -0.3450 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2815 -0.2806 -0.2827 -0.2800 -0.2796
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2723 -0.2714 -0.2718 -0.2708 -0.2703
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3493 -0.3488 -0.3483 -0.3470 -0.3470 -0.3111
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110  0.06506 -0.05055  0.07451  0.07939
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3344 -0.3341 -0.3338 -0.3335 -0.3333 -0.3111
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110  0.05491 -0.06203  0.06138  0.06807
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3111 -0.3111 -0.3111 -0.3111 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2853 -0.2849 -0.2874 -0.2844 -0.2841
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2661 -0.2651 -0.2655 -0.2644 -0.2638
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3107 -0.3107 -0.3108 -0.3107 -0.3107
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110  0.29952  0.03818  0.30757  0.31579
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3227 -0.3226 -0.3224 -0.3223 -0.3222 -0.3111
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110  0.02354 -0.07831  0.03184  0.03613
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4206 -0.4190 -0.4173 -0.4165 -0.4149 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2849 -0.2843 -0.2846 -0.2839 -0.2835
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2918 -0.2914 -0.2915 -0.2911 -0.2908
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2433 -0.2419 -0.2424 -0.2409 -0.2399
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3230 -0.3229 -0.3227 -0.3226 -0.3224 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4168 -0.4153 -0.4137 -0.4128 -0.4114 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2942 -0.2939 -0.2940 -0.2936 -0.2934
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2917 -0.2912 -0.2914 -0.2909 -0.2907
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2800 -0.2794 -0.2797 -0.2789 -0.2785
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110 -0.31110  0.05324  0.47268  0.48811
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2967 -0.2963 -0.2965 -0.2961 -0.2959
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4446 -0.4426 -0.4405 -0.4384 -0.4377 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3143 -0.3143 -0.3142 -0.3142 -0.3142 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3335 -0.3331 -0.3328 -0.3324 -0.3322 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2496 -0.2479 -0.2503 -0.2469 -0.2460
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3716 -0.3708 -0.3698 -0.3680 -0.3686 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.4320 -0.4302 -0.4284 -0.4267 -0.4258 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3009 -0.3006 -0.3007 -0.3005 -0.3003
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    # -0.31110 -0.31110 -0.04972 -0.11272 -0.04163 -0.03742
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3154 -0.3154 -0.3153 -0.3153 -0.3152 -0.3111
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.2862 -0.2852 -0.2867 -0.2849 -0.2845
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # -0.3111 -0.3095 -0.3095 -0.3095 -0.3094 -0.3094
}

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

## Compute some reduced dims
set.seed(20220423)
sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 20)
sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

## We don't want to model the pathology groups as integers / numeric
## so let's double check this
stopifnot(is.factor(sce_pseudo$path_groups))

## Add APOe genotype info
sce_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E3", "Br3880" = "E3/E3", "Br3874" = "E2/E3")[sce_pseudo$subject]

## For the spatialLIBD shiny app
rowData(sce_pseudo)$gene_search <-
    paste0(
        rowData(sce_pseudo)$gene_name,
        "; ",
        rowData(sce_pseudo)$gene_id
    )

## Drop things we don't need
spatialCoords(sce_pseudo) <- NULL
imgData(sce_pseudo) <- NULL

## Simplify the colData()  for the pseudo-bulked data
colData(sce_pseudo) <- colData(sce_pseudo)[, sort(c(
    "age",
    "sample_id",
    "path_groups",
    "subject",
    "sex",
    "pmi",
    "APOe",
    "race",
    "diagnosis",
    "rin",
    "BCrating",
    "braak",
    "cerad"
))]

## Load pathology colors
## This info is used by spatialLIBD v1.7.18 or newer
source(here("code", "colors_pathology.R"), echo = TRUE, max.deparse.length = 500)
sce_pseudo$path_groups_colors <- colors_pathology[as.character(sce_pseudo$path_groups)]

## save RDS file
saveRDS(
    sce_pseudo,
    file = file.path(
        dir_rdata,
        paste0("sce_pseudo_pathology_", opt$spetype, ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
