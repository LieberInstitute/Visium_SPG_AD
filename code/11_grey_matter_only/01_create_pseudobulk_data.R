# library(sgejobs)

# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "01_create_pseudobulk_data",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "10G")
#To execute the script builder, use: sh 01_create_pseudobulk_data.sh

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



library("SpatialExperiment")
library("here")
library("spatialLIBD")
library("scuttle")
library("edgeR")
library("sessioninfo")

##output directory
dir_rdata<- here::here("processed-data","11_grey_matter_only", opt$spetype)
dir.create(dir_r, showWarnings = FALSE)

##load spe data
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            opt$spetype,
            paste0("spe_harmony_",opt$spetype, ".rds")

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

spe <-cluster_import(
    spe, cluster_dir = here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
),
prefix = "")

##subset spe data based on subject and cluster 1 for k = 2
spe_new <- spe[,!spe$subject == "Br3874"]

if(opt$spetype == "wholegenome") {
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


##pseudobulk across pathology labels
sce_pseudo <- aggregateAcrossCells(
    spe_new,
    DataFrame(path_groups = spe_new$path_groups,
              sample_id = spe_new$sample_id
    )
)
x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)

stopifnot(identical(rownames(x), rownames(sce_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(sce_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <- x

##save RDS file
saveRDS(
    sce_pseudo,
    file = here::here(
        "processed-data","11_grey_matter_only", opt$spetype,
        paste0("sce_pseudo_path_type",opt$spetype, ".RDS")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




