# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "01_spatial_registration",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "10G",
#     command = "Rscript 01_spatial_registration.R",
#     create_logdir = TRUE
#
# )


library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)
library(sessioninfo)

# ##create directories
# dir.create(here::here("processed-data","10_spatial_registration", spetype), showWarnings = FALSE)
# dir.create(here::here("processed-data","10_spatial_registration", "pseudo_bulked", spetype), showWarnings = FALSE)
# dir.create(here::here("processed-data", "10_spatial_registration", "dupCor", spetype), showWarnings = FALSE)
# dir.create(here::here("processed-data","10_spatial_registration", "specific_Ts", spetype), showWarnings = FALSE)
# dir.create(here::here("code", "10_spatial_registration"))

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)
spetype <- commandArgs(trailingOnly = TRUE)
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            spetype,
            paste0("spe_harmony_", spetype, ".rds")
        )
    )

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        spetype,
        "clusters_BayesSpace"
    ),
    prefix = ""
)

sce_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        BayesSpace = colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]],
        sample_id = spe$sample_id
    )
)

x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)

stopifnot(identical(rownames(x), rownames(sce_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(sce_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)

rowData(sce_pseudo)$low_expr <- filterByExpr(sce_pseudo)
summary(rowData(sce_pseudo)$low_expr)

# saveRDS(
#     sce_pseudo,
#     file = here::here(
#         "processed-data",
#         "10_spatial_registration",
#         "pseudo_bulked",
#         paste0("sce_pseudobulked_BayesSpace", k_nice, spetype, ".RDS")
#     )
# )

###############################
##### get mean expression  ####
mat <- logcounts(sce_pseudo)

#####################
## Build a group model

# convert variables to factors
sce_pseudo$spatial.cluster <-
    as.factor(colData(sce_pseudo)[[paste0("BayesSpace_harmony_k", k_nice)]])
table(is.na(sce_pseudo$spatial.cluster))
sce_pseudo$sex <- as.factor(sce_pseudo$sex)
sce_pseudo$diagnosis <- as.factor(sce_pseudo$diagnosis)
sce_pseudo$sample_id <- as.factor(sce_pseudo$sample_id)
sce_pseudo$age <- as.numeric(sce_pseudo$age)


mod <- with(
    colData(sce_pseudo),
    model.matrix(~ 0 + spatial.cluster + diagnosis + age + sex)
)

colnames(mod) <- gsub("cluster", "", colnames(mod))

## get duplicate correlation
# the design matrix of the micro-array experiment, with rows
# corresponding to arrays and columns to comparisons to be
# estimated. The number of rows must match the number of
# columns of ‘object’.

corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudo$sample_id
)
message("Detected correlation: ", corfit$consensus.correlation)

# > dim(mod)
# [1] 21  7
# > dim(mat)
# [1] 18116   100

# saveRDS(
#     corfit,
#     file = here::here(
#         "processed-data",
#         "10_spatial_registration",
#         "dupCor",
#         paste0("pseudobulked_dupCor_k", k_nice, spetype ,".RDS")
#     )
# )

## Next for each layer test that layer vs the rest
cluster_idx <-
    splitit(sce_pseudo$spatial.cluster)

eb0_list_cluster <- lapply(cluster_idx, function(x) {
    res <- rep(0, ncol(sce_pseudo))
    res[x] <-
        1 # indicator of whether pseudobulked column belongs to
    m <-
        with(
            colData(sce_pseudo),
            # find genes that are diff expressed across diff BayesSpace clusters adjusting for
            model.matrix(~ res + diagnosis + age + sex)
        ) # age, diagnosis, sex
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudo$sample_id,
            ## block helps take into account that there are some bulks that come from the same sample and
            ## they may be highly correlated to each other
            correlation = corfit$consensus.correlation
        )
    )
})


# saveRDS(
#     eb0_list_cluster,
#     file = here::here(
#         "processed-data",
#         "10_spatial_registration",
#         "specific_Ts",
#         paste0("pseudobulked_specific_Ts_k", k_nice, spetype,".RDS")
#     )
# )
#

##########
## Extract the p-values

pvals0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cluster) <- rownames(mat)

t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cluster) <- rownames(mat)
fdrs0_contrasts_cluster <- apply(pvals0_contrasts_cluster, 2, p.adjust, "fdr")

data.frame(
    "FDRsig" = colSums(fdrs0_contrasts_cluster < 0.05 &
        t0_contrasts_cluster > 0),
    "Pval10-6sig" = colSums(pvals0_contrasts_cluster < 1e-6 &
        t0_contrasts_cluster > 0),
    "Pval10-8sig" = colSums(pvals0_contrasts_cluster < 1e-8 &
        t0_contrasts_cluster > 0)
)

# for k = 04
# FDRsig Pval10.6sig Pval10.8sig
# 1  16576        2113         401
# 2    266          12           1
# 3     38           0           0
# 4      0           0           0

###################

ground_truth <- spatialLIBD::fetch_data("modeling_results")

cor_stats_layer <- layer_stat_cor(
    t0_contrasts_cluster,
    modeling_results = ground_truth,
    model_type = "enrichment",
    top_n = 100
)

## plot output directory
dir_plots <-
    here::here("plots", "10_spatial_registration", spetype)
# dir.create(dir_plots, showWarnings = FALSE)

# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html newer function for plotting

pdf(
    file = here::here(
        "plots",
        "10_spatial_registration",
        paste0(
            "enrichment_analysis_k_",
            k_nice, "_",
            spetype,
            ".pdf"
        )
    ),
    width = 8
)
layer_stat_cor_plot(
    cor_stats_layer,
    max = 1
)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
