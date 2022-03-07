#library(remotes)
#remotes::install_github('LieberInstitute/sgejobs')
# library(sgejobs)
# sgejobs::job_loop(
#     loops = list(spetype = c(
#         "wholegenome", "targeted"
#     )),
#     name = "fasthplus_optimal_k",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     task_num = 15,
#     tc = 14
#
# )


## Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "spetype", "s", 2, "character", "SPE type: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}


library('here')
#remotes::install_github(repo="ntdyjack/fasthplus")
library('fasthplus')
library("sessioninfo")
library("spatialLIBD")

if (FALSE) {
    ## For testing, run interactively if needed
    opt <- list(spetype = "wholegenome")
}

## Load the data
spe <- readRDS(here::here(
    "processed-data",
    "08_harmony_BayesSpace",
    opt$spetype,
    paste0("spe_harmony_", opt$spetype, ".rds")
))

# import cluster info
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

##create output directories
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)


##choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)


##Code from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/06_fastplus/06_fasthplus.R
#hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
set.seed(20220304)
## Use 5% of the data
fasthplus <- hpb(D= reducedDims(spe)$HARMONY,L=colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]],t=floor(nrow(spe) * 0.05),r=30)
results <- data.frame (k=k, fasthplus=fasthplus, type = opt$spetype)
write.table(results, file = here::here("processed-data","08_harmony_BayesSpace","fasthplus_results.csv"), append = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


