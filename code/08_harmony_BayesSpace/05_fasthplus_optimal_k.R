# library(remotes)
# remotes::install_github('LieberInstitute/sgejobs')
# library(sgejobs)
# sgejobs::job_loop(
#     loops = list(spetype = c(
#         "wholegenome", "targeted"
#     )),
#     name = "fasthplus_optimal_k",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     task_num = 28
# )


## Required libraries
library("getopt")

## Specify parameters
spec <- matrix(
    c(
        "spetype",
        "s",
        2,
        "character",
        "SPE type: wholegenome or targeted",
        "help",
        "h",
        0,
        "logical",
        "Display help"
    ),
    byrow = TRUE,
    ncol = 5
)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}


library("here")
# remotes::install_github(repo="ntdyjack/fasthplus")
library("fasthplus")
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

## create output directories
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)


## choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

## Find the right value of t to use
## Adapted code from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/06_fastplus/06_fasthplus.R#L46-L63
find_t <- function(L, proportion = 0.05) {
    t_value <- floor(length(L) * proportion)
	message(Sys.time(), " t value based on proportion: ", t_value)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    t_value <- ifelse(smallest_cluster_size > (t_value / n_labels), t_value, smallest_cluster_size * n_labels)
	message(Sys.time(), " t value estimated: ", t_value)
	return(t_value)
}

## Function that takes a logical vector of spots we can use
## and returns the list of those we can actually use and the resulting t
find_spots_and_t <- function(usable_spots, t_proportion = 0.05) {
	L <- colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]]
	t_value <- find_t(L = L[usable_spots], proportion = t_proportion)

	cluster_prop <- table(L[usable_spots]) / sum(usable_spots)
	bad_clusters <- which(cluster_prop < t_proportion / k)
	if(length(bad_clusters) > 0) {
	    message("For k: ", k, " we are dropping small clusters: ", paste(names(bad_clusters), collapse = ", "))
		usable_spots[ L %in% as.integer(names(bad_clusters)) ] <- FALSE
	    t_value <- find_t(L = L[usable_spots], proportion = t_proportion)
	}
	return(list(t_value = t_value, spots_to_use = usable_spots))
}
spots_and_t <- find_spots_and_t(usable_spots = rep(TRUE, ncol(spe)))


## Code from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/06_fastplus/06_fasthplus.R
# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
set.seed(20220304)
fasthplus <-
    hpb(
        D = reducedDims(spe[, spots_and_t$spots_to_use])$HARMONY,
        L = colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]][spots_and_t$spots_to_use],
        t = spots_and_t$t_value,
        r = 30
    )
results <-
    data.frame(
        k = k,
        fasthplus = fasthplus,
        type = opt$spetype,
		spots_set = "all_spots",
		t_value = spots_and_t$t_value
    )
write.table(
    results,
    file = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "fasthplus_results.csv"
    ),
    append = TRUE,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

## Repeat but after dropping the white matter
if(opt$spetype == "wholegenome") {
	spots_and_t_noWM <- find_spots_and_t(usable_spots = spe$BayesSpace_harmony_k02 != 2)
} else {
	spots_and_t_noWM <- find_spots_and_t(usable_spots = spe$BayesSpace_harmony_k04 != 4)
}


set.seed(20220304)
fasthplus <-
    hpb(
        D = reducedDims(spe[, spots_and_t_noWM$spots_to_use])$HARMONY,
        L = colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]][spots_and_t_noWM$spots_to_use],
        t = spots_and_t_noWM$t_value,
        r = 30
    )
results <-
    data.frame(
        k = k,
        fasthplus = fasthplus,
        type = opt$spetype,
		spots_set = "grey_matter",
		t_value = spots_and_t_noWM$t_value
    )
write.table(
    results,
    file = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "fasthplus_results.csv"
    ),
    append = TRUE,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
