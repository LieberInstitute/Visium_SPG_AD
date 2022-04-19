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

## Extract the data
mat <- assays(sce_pseudo)$logcounts

## Compute correlation
## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/7973fcebb7c4b17cc3e23be2c31ac324d1cc099b/code/10_spatial_registration/01_spatial_registration.R#L134-L150
mod <- with(colData(sce_pseudo),
    model.matrix( ~ 0 + path_groups + age + sex))

corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudo$sample_id)
message("Detected correlation: ", corfit$consensus.correlation)

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
patho_idx <- splitit(sce_pseudo$path_groups)

eb0_list <- lapply(patho_idx, function(x) {
    res <- rep(0, ncol(sce_pseudo))
    res[x] <- 1
    m <- model.matrix(~ res + age + sex)
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudo$sample_id,
            correlation = corfit$consensus.correlation
        )
    )
})

######### PAIRWISE t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1355-L1383

## Build a group model
mod <- with(colData(sce_pseudo), model.matrix(~ 0 + path_groups))

fit <-
    lmFit(
        mat,
        design = mod,
        block = sce_pseudo$sample_id,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)


## Define the contrasts for each layer vs the rest (excluding WM comparisons since we have that one already)
path_combs <- combn(colnames(mod), 2)
path_contrasts <- apply(path_combs, 2, function(x) {
    z <- paste(x, collapse = '-')
    makeContrasts(contrasts = z, levels = mod)
})
rownames(path_contrasts) <- colnames(mod)
colnames(path_contrasts) <-
    apply(path_combs, 2, paste, collapse = '-')
eb_contrasts <- eBayes(contrasts.fit(fit, path_contrasts))


######### ANOVA t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity_fstats.R#L24-L85

## From layer_specificity.R
fit_f_model <- function(sce) {
    message(paste(Sys.time(), 'starting the model run'))

    ## Extract the data
    mat <- assays(sce)$logcounts

    ## For dropping un-used levels
    sce$path_groups <- factor(sce$path_groups)

    ## Build a group model
    mod <- with(colData(sce), model.matrix( ~ path_groups + age + sex))

    ## Takes like 2 min to run
    corfit <-
        duplicateCorrelation(mat, mod, block = sce$subject)
    message(paste(Sys.time(), 'correlation:', corfit$consensus.correlation))
    fit <-
        lmFit(
            mat,
            design = mod,
            block = sce$subject,
            correlation = corfit$consensus.correlation
        )
    eb <- eBayes(fit)
    return(eb)
}

ebF_list <-
    lapply(list('noWM' = sce_pseudo), fit_f_model)

## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
    x <- ebF_list[[i]]
    top <-
        topTable(
            x,
            coef = grep("path_groups", colnames(x$coefficients)),
            sort.by = 'none',
            number = length(x$F)
        )
    # identical(p.adjust(top$P.Value, 'fdr'), top$adj.P.Val)
    res <- data.frame(
        'f' = top$F,
        'p_value' = top$P.Value,
        'fdr' = top$adj.P.Val,
        'AveExpr' = top$AveExpr,
        stringsAsFactors = FALSE
    )
    colnames(res) <- paste0(i, '_', colnames(res))
    return(res)
}))
f_stats$ensembl <- rownames(sce_pseudo)
f_stats$gene <- rowData(sce_pseudo)$gene_name
rownames(f_stats) <- NULL

head(f_stats)


save(
    f_stats,
    eb0_list,
    eb_contrasts,
    file = file.path(dir_rdata, "pathology_model_results_raw.Rdata")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
