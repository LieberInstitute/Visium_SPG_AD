library("here")
library("sessioninfo")
library("spatialLIBD")
library("EnhancedVolcano")

## Locate data directory
dir_rdata <- here::here(
        "code",
        "05_deploy_app_wholegenome"
    )

load(file.path(dir_rdata, "Visium_IF_AD_modeling_results.Rdata"),
    verbose = TRUE
)
sce_pseudo <-
    readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))


## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

## Build a data.frame with the info needed first
## Adapted from https://github.com/LieberInstitute/spatial_hpc/blob/62334e88788c2a7010f2d500418639769121ecf1/code/08_pseudobulk/PRECAST/plot_volcano.R#L75-L80

testname <- "Ab"
make_volcano <- function(testname) {
    subset_sig <- subset(sig_genes, test == testname & model_type == "enrichment")
    df <- data.frame(
        gene = subset_sig$gene,
        logFC = subset_sig$logFC,
        FDR = subset_sig$fdr,
        sig = subset_sig$fdr < 0.05
    )

    ## Adapted from https://github.com/LieberInstitute/spatial_hpc/blob/62334e88788c2a7010f2d500418639769121ecf1/code/08_pseudobulk/PRECAST/plot_volcano.R#L580-L591
    EnhancedVolcano(df,
        lab = df$gene,
        x = 'logFC',
        y = 'FDR',
        FCcutoff = 1,
        pCutoff = 0.05,
        ylab = "-log10 FDR",
        legendLabels = c('Not sig.','Log (base 2) FC','FDR',
            'FDR & Log (base 2) FC'),
        title = paste(testname, "vs others"),
        subtitle = ""
    )
}



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
