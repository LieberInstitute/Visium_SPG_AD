# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "03_model_pathology",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 03_model_pathology.sh

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


library("here")
library("sessioninfo")
library("SpatialExperiment")
library("rafalib")
library("limma")


## output directory
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here::here("plots", "11_grey_matter_only", opt$spetype)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load spe data
sce_pseudo <-
    readRDS(
        here::here(
            "processed-data",
            "11_grey_matter_only",
            opt$spetype,
            paste0("sce_pseudo_pathology_", opt$spetype, ".rds")
        )
    )

## We don't want to model the pathology groups as integers / numeric
## so let's double check this
stopifnot(is.factor(sce_pseudo$path_groups))

## Add APOe genotype info
sce_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E4", "Br3800" = "E3/E3", "Br3874" = "E2/E3")[sce_pseudo$subject]

## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
layer_idx <- splitit(sce_layer$layer_guess)

eb0_list <- lapply(layer_idx, function(x) {
    res <- rep(0, ncol(sce_layer))
    res[x] <- 1
    m <- model.matrix(~ res)
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_layer$subject_position,
            correlation = corfit$consensus.correlation
        )
    )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
