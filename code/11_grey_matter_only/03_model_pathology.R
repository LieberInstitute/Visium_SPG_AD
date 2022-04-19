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
library("spatialLIBD")
library("scater")
library("jaffelab")


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

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68
pca <- prcomp(t(assays(sce_pseudo)$logcounts))
message(Sys.time(), " % of variance explained for the top 50 PCs:")
jaffelab::getPcaVars(pca)[seq_len(50)]
pca_pseudo<- pca$x[, seq_len(50)]
reducedDims(sce_pseudo) <- list(PCA=pca_pseudo)

## Plot PCs with different colors
## Each point here is a sample
pdf(file = file.path(dir_plots, paste0("sce_pseudo_pca.pdf")), width = 14, height = 14)
plotPCA(sce_pseudo, colour_by = "age", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "path_groups", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "subject", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "pmi", ncomponents = 12, point_size = 1)
plotPCA(sce_pseudo, colour_by = "APOe", ncomponents = 12, point_size = 1)
dev.off()


## Obtain percent of variance explained at the gene level
## using scater::getVarianceExplained()
vars <- getVarianceExplained(sce_pseudo,
    variables = c(
        "age",
        "sample_id",
        "path_groups",
        "subject",
        "sex",
        "pmi",
        "APOe"
    )
)

## Now visualize the percent of variance explained across all genes
pdf(file = file.path(dir_plots, paste0("sce_pseudo_gene_explanatory_vars.pdf")))
plotExplanatoryVariables(vars)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
