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

# Ensure that necessary packages are installed
packages_needed <- c("jaffelab", "devtools", "remotes") # add other packages as needed
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]

if (length(packages_to_install) > 0) {
    # Install remotes if not already installed
    if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
    }
    
    # Install missing packages from CRAN
    cran_packages <- packages_to_install[packages_to_install %in% rownames(available.packages())]
    if (length(cran_packages) > 0) {
        install.packages(cran_packages)
    }
    
    # Install jaffelab from GitHub
    if ("jaffelab" %in% packages_to_install) {
        remotes::install_github("LieberInstitute/jaffelab")
    }
    
    # Load all required packages
    invisible(lapply(packages_needed, library, character.only = TRUE))
} else {
    message("All packages are already installed.")
}

# Check and install missing packages
necessary_packages <- c("getopt", "here", "spatialLIBD", "sessioninfo", "scater")
missing_packages <- necessary_packages[!(necessary_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

# Load packages
lapply(necessary_packages, library, character.only = TRUE)

# Required libraries
library("getopt")

## Specify parameters
#spec <- matrix(c(
#    "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
#    "help", "h", 0, "logical", "Display help"
#), byrow = TRUE, ncol = 5)
#opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
#if (!is.null(opt$help)) {
#    cat(getopt(spec, usage = TRUE))
#    q(status = 1)
#}

# Manually define the options instead of using command line arguments
opt <- list(spetype = "wholegenome")  # Change "wholegenome" to "targeted" if needed

# Now proceed with using opt$spetype to set directory paths
library("here")  # Ensure 'here' library is loaded for path management
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
print(dir_rdata)  # This will print the directory path to confirm it's correctly set

# Create the directory if it doesn't exist
if (!dir.exists(dir_rdata)) {
    dir.create(dir_rdata, recursive = TRUE)
    cat("Directory created: ", dir_rdata, "\n")
} else {
    cat("Directory already exists: ", dir_rdata, "\n")
}

## For testing
if (FALSE) {
    opt <- list(spetype = "wholegenome")
}

library("here")
library("spatialLIBD")
stopifnot(packageVersion("spatialLIBD") >= "1.9.18")
library("sessioninfo")
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
            "Ab",
            "n_Ab",
            "pTau",
            "n_pTau",
            "both",
            "n_both"
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

summary(spe$path_groups)
stopifnot(is.factor(spe$path_groups))


## pseudobulk across pathology labels
sce_pseudo <-
    registration_pseudobulk(spe,
        var_registration = "path_groups",
        var_sample_id = "sample_id",
        min_ncells = 15
    )
dim(sce_pseudo)

## Add APOe genotype info
sce_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E3", "Br3880" = "E3/E3", "Br3874" = "E2/E3")[sce_pseudo$subject]

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
    "cerad",
    "ncells"
))]

## Explore the resulting data
options(width = 400)
as.data.frame(colData(sce_pseudo))

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

## For the spatialLIBD shiny app
rowData(sce_pseudo)$gene_search <-
    paste0(
        rowData(sce_pseudo)$gene_name,
        "; ",
        rowData(sce_pseudo)$gene_id
    )

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
