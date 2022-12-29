# library("sgejobs")

# sgejobs::job_single(
#     "mathys",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 03_mathys.R"
# )

#### load relevant packages ####
library("readxl")
library("spatialLIBD")
library("dplyr")
library("sessioninfo")
library("here")
library("scran")
library("purrr")

### load get_ensemble function
source(here("code/14_external_gene_sets/get_ensembl_function.R"))

### load modeling results
load(here(
    "processed-data", "11_grey_matter_only", "wholegenome",
    "Visium_IF_AD_modeling_results.Rdata"
))

# Number of sets: 6 cell types * 4 models = 24 sets * direction (2) = 48 sets
#
# Note: snRNA-seq
#
# Direction available: IndModel.FC or MixedModel.z for the indicator model (Wilcoxon test) or for the mixed effects model (Poisson)
# Statistics available: Yes, but they also took into account effect size with DEGs.Ind.Model and DEGs.Ind.Mix.models
# Ignore direction: 6 cell types * 4 models = 24 sets
# Use direction: 6 cell types * 4 models * 2 directions = 48 sets
# For each of the above, check how many genes we have in each set. They might be too small.
# Note that Excel might introduce some issues in this table too.


# no-pathology vs pathology differential expression = np vs p, 1
# no-pathology vs early-pathology differential expression = np vs ep, 12
# early-pathology vs late-pathology differential expression = ep vs lp, 23

input_dir <- here("raw-data", "GeneSets", "2_snRNA-seq",
            "1_Mathys et al_PFC", "Mathys et al.xlsx")

mathys_ex <- read_excel(input_dir, sheet = "Ex", skip = 1)
mathys_in <- read_excel(input_dir, sheet = "In", skip = 1)
mathys_ast <- read_excel(input_dir, sheet = "Ast", skip = 1)
mathys_oli <- read_excel(input_dir, sheet = "Oli", skip = 1)
mathys_opc <- read_excel(input_dir, sheet = "Opc", skip = 1)
mathys_mic <- read_excel(input_dir, sheet = "Mic", skip = 1)


change_col_names <- function(df) {
    colnames(df)[1] <- "gene_set_np_v_p"
    colnames(df)[12] <- "gene_set_np_v_ep"
    colnames(df)[23] <- "gene_set_ep_v_lp"
    df
}


df.list <- list(mathys_ex, mathys_in, mathys_ast, mathys_oli, mathys_opc, mathys_mic)
df.list <- purrr::map(df.list, change_col_names)


res_1 <- purrr::map(df.list, get_ensembl, gene_set_np_v_p, "gene_set_np_v_p")
res_2 <- purrr::map(df.list, get_ensembl, gene_set_np_v_ep, "gene_set_np_v_ep")
res_3 <- purrr::map(df.list, get_ensembl, gene_set_ep_v_lp, "gene_set_ep_v_lp")

# mathys_ex, mathys_in, mathys_ast, mathys_oli, mathys_opc, mathys_mic

mathys_geneList <- list(
    np_v_p_ex = res_1[[1]]$gene_ensembl_id,
    np_v_p_in = res_1[[2]]$gene_ensembl_id,
    np_v_p_ast = res_1[[3]]$gene_ensembl_id,
    np_v_p_oli = res_1[[4]]$gene_ensembl_id,
    np_v_p_opc = res_1[[5]]$gene_ensembl_id,
    np_v_p_mic = res_1[[6]]$gene_ensembl_id,
    np_v_ep_ex = res_2[[1]]$gene_ensembl_id,
    np_v_ep_in = res_2[[2]]$gene_ensembl_id,
    np_v_ep_ast = res_2[[3]]$gene_ensembl_id,
    np_v_ep_oli = res_2[[4]]$gene_ensembl_id,
    np_v_ep_opc = res_2[[5]]$gene_ensembl_id,
    np_v_ep_mic = res_2[[6]]$gene_ensembl_id,
    ep_v_lp_ex = res_3[[1]]$gene_ensembl_id,
    ep_v_lp_in = res_3[[2]]$gene_ensembl_id,
    ep_v_lp_ast = res_3[[3]]$gene_ensembl_id,
    ep_v_lp_oli = res_3[[4]]$gene_ensembl_id,
    ep_v_lp_opc = res_3[[5]]$gene_ensembl_id,
    ep_v_lp_mic = res_3[[6]]$gene_ensembl_id
)

mathys_enrichment <- gene_set_enrichment(
    mathys_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

mathys_depleted <- gene_set_enrichment(
    mathys_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)

##### enriched genes plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/03_mathys_enriched.pdf"), width = 15)
gene_set_enrichment_plot(
    mathys_enrichment,
    xlabs = unique(mathys_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(mathys_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

##### depleted genes plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/03_mathys_depleted.pdf"), width = 15)
gene_set_enrichment_plot(
    mathys_depleted,
    xlabs = unique(mathys_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(mathys_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
