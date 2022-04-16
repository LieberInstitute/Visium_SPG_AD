library(sgejobs)
sgejobs::job_loop(
   loops = list(speopt$spetype = c(
       "wholegenome", "targeted"
    )),
    name = "01_create_pseudobulk_data",
    create_shell = TRUE,
    queue = "bluejay",
    memory = "10G")
#To execute the script builder, use: sh 01_create_pseudobulk_data.sh



# Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "speopt$spetype", "s", 2, "character", "SPE opt$spetype: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}



library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)



dir.create(here::here("processed-data","gm_only", opt$spetype), showWarnings = FALSE)

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


spe_new <- spe[,!spe$subject == "Br3874"]
spe_new <- spe_new[, !spe_new$BayesSpace_harmony_k02 == 2]
# > dim(colData(spe))
# [1] 38115   111
# # > dim(colData(spe_new))
# [1] 21086   111

## > unique(spe$path_groups)
# [1] "none"      "next_Ab+"  "Ab+"       "both"      "pT+"       "next_both"
# [7] "next_pT+"


    ##pseudobulk across pathology labels
    ##.Rdata stored here for path levels
    #     here::here(
    #         "processed-data",
    #         "09_pathology_vs_BayesSpace",
    #         "pathology_levels")



