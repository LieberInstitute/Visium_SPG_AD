
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

# load("sce_for_iSEE_LS.rda", verbose = TRUE)# load the pseudobulked object sce_IF
sce_IF <- readRDS("sce_pseudo_pathology_wholegenome.rds")

## Make unique gene names
rownames(sce_IF) <-
    uniquifyFeatureNames(rowData(sce_IF)$gene_id, rowData(sce_IF)$gene_name)

# stopifnot(all(unique(sce_IF$BayesSpace) %in% names(cell_cols.clean)))

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce_IF)
# 876.33 MB

source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce_IF) <- cbind(
    colData(sce_IF)[, !colnames(colData(sce_IF)) %in% c("subject", "path_groups")],
    colData(sce_IF)[, c("path_groups", "subject")]
)

sce_IF$subject <- as.factor(sce_IF$subject)

sce_IF <- registerAppOptions(sce_IF, color.maxlevels = length(colData(sce_IF)$path_groups_colors))

iSEE(
    sce_IF,
    appTitle = "Kwon2022_pseudobulk_AD_pathology_wholegenome",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
    subject = function(n) {
        cols <- paletteer::paletteer_d(
            palette = "RColorBrewer::Dark2",
            n = length(unique(sce_IF$subject))
        )
        cols <- as.vector(cols)
        names(cols) <- levels(sce_IF$subject)
        return(cols)
    },
        path_groups = function(n) {
            return(colData(sce_IF)$path_groups_colors)
        }
    ))
)
