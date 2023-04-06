

# Number of sets: 6 cell types * concordance (2 options) = 12 sets + direction 6 sets = 18 sets
#
# Note: snRNA-seq
#
# Direction available: Grubman.LogFC (AD vs Control a priori). There's not that many genes.
#
# Statistics available: No. But we have Concordance.
## Correction: We have concordance only for the Mathys comparison. FDR available all others.
#
#  For each cell type, use all genes (ignore Concordance)
#  For each cell type, filter to Concordance == TRUE

### not enough genes with Concordance == FALSE


#  For each cell type, use all genes (ignore Concordance), and separate by direction.

# library("sgejobs")

# sgejobs::job_single(
#     "grubman",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 04_grubman.R"
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


input_dir <- here("raw-data", "GeneSets", "2_snRNA-seq",
                  "2_Grubman et al_Entorhinal", "Grubman et al.xlsx")

# comparison with Mathys
table_s3_grubman <- read_excel(input_dir,
    sheet = "Supplementary Table 3", skip = 10,
    col_names = TRUE
)

##renaming column for readability
table_s3_grubman$grubman_logfc <- table_s3_grubman$`Grubman.LogFC (AD vs Control a priori)`
table_s3_grubman <- table_s3_grubman |> subset(select = -`Grubman.LogFC (AD vs Control a priori)`)



#### NOTES ON TABLES ####

# > nrow(table_s3_grub)
# [1] 94                                                     col_names = TRUE)

# > colnames(table_s3_grub)
# [1] "...1"                                    "Genes"
# [3] "Concordance"                             "Grubman.LogFC (AD vs Control a priori)"
# [5] "Mathys.LogFC(no pathology vs pathology)" "cell type"

# > unique(table_s3_grub$`cell type`)
# [1] "astro"               "mg"                  "neuron (excitatory)"
# [4] "neuron (inhibitory)" "oligo"               "OPC"





#### ALL GENES  ####
table_s3_grubman_astro <- table_s3_grubman |> dplyr::filter(`cell type` == "astro")
table_s3_grubman_mg <- table_s3_grubman |> dplyr::filter(`cell type` == "mg")
table_s3_grubman_ex <- table_s3_grubman |> dplyr::filter(`cell type` == "neuron (excitatory)")
table_s3_grubman_inh <- table_s3_grubman |> dplyr::filter(`cell type` == "neuron (inhibitory)")
table_s3_grubman_oligo <- table_s3_grubman |> dplyr::filter(`cell type` == "oligo")
table_s3_grubman_OPC <- table_s3_grubman |> dplyr::filter(`cell type` == "OPC")

# > nrow(table_s3_grubman_astro )
# [1] 32
# > nrow(table_s3_grubman_mg)
# [1] 11
# > nrow(table_s3_grubman_ex)
# [1] 21
# > nrow(table_s3_grubman_inh )
# [1] 4
# > nrow(table_s3_grubman_oligo)
# [1] 22
# > nrow(table_s3_grubman_OPC)
# [1] 4



df_mathys_list <- list(
    table_s3_grubman_astro, table_s3_grubman_mg, table_s3_grubman_ex, table_s3_grubman_inh,
    table_s3_grubman_oligo, table_s3_grubman_OPC
)

res_1 <- purrr::map(df_mathys_list, get_ensembl, Genes, "Genes")


grubman_geneList <- list(
    grubman_mathys_astro = res_1[[1]],
    grubman_mathys_mg = res_1[[2]],
    grubman_mathys_ex = res_1[[3]],
    grubman_mathys_inh = res_1[[4]],
    grubman_mathys_oligo = res_1[[5]],
    grubman_mathys_OPC = res_1[[6]]
)



grubman_enrichment <- gene_set_enrichment(
    grubman_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

dummy_row_enrichment <- data.frame(OR = 2,
                                   Pval = 0.01 ,
                                   test = "none" ,
                                   NumSig = 0  ,
                                   SetSize = 0 ,
                                   ID = "dummy",
                                   model_type = "enrichment" ,
                                   fdr_cut = 0.1)

grubman_enrichment  = rbind(grubman_enrichment,
                              dummy_row_enrichment)

print("grubman_enrichment")
print(grubman_enrichment)

grubman_depleted <- gene_set_enrichment(
    grubman_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)

dummy_row_depleted <- data.frame(OR = 2,
                                 Pval = 0.01 ,
                                 test = "none" ,
                                 NumSig = 0  ,
                                 SetSize = 0 ,
                                 ID = "dummy",
                                 model_type = "depletion" ,
                                 fdr_cut = 0.1)
grubman_depleted  = rbind(grubman_depleted,
                            dummy_row_depleted)


print("grubman_depleted")
print(grubman_depleted)

##### enrichment plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/04_grubman_enriched.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_enrichment,
    xlabs = unique(grubman_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()


pdf(paste0(output_dir, "/04_grubman_depleted.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_depleted,
    xlabs = unique(grubman_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()


# #### GENES SPLIT BY DIRECTION  ####
#
# table_s3_grubman_astro_pos <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "astro" & grubman_logfc > 0)
#
# table_s3_grubman_mg_pos <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "mg"& grubman_logfc > 0)
#
# table_s3_grubman_ex_pos  <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "neuron (excitatory)" & grubman_logfc > 0)
#
# table_s3_grubman_inh_pos <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "neuron (inhibitory)" & grubman_logfc > 0)
#
# table_s3_grubman_oligo_pos <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "oligo" & grubman_logfc > 0)
#
# table_s3_grubman_OPC_pos <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "OPC" & grubman_logfc > 0)
#
# nrow(table_s3_grubman_astro_neg )
#
# nrow(table_s3_grubman_mg_neg)
#
# nrow(table_s3_grubman_ex_neg)
#
# nrow(table_s3_grubman_inh_neg )
#
# nrow(table_s3_grubman_oligo_neg)
#
# nrow(table_s3_grubman_OPC_neg)
#
# # > nrow(table_s3_grubman_astro_pos )
# # [1] 17
# # > nrow(table_s3_grubman_mg_pos)
# # [1] 9
# # > nrow(table_s3_grubman_ex_pos)
# # [1] 14
# # > nrow(table_s3_grubman_inh_pos )
# # [1] 2
# # > nrow(table_s3_grubman_oligo_pos)
# # [1] 18
# # > nrow(table_s3_grubman_OPC_pos)
# # [1] 3
#
# # > nrow(table_s3_grubman_astro_neg )
# # >
# #     [1] 15
# # > nrow(table_s3_grubman_mg_neg)
# # [1] 2
# # >
# #     > nrow(table_s3_grubman_ex_neg)
# # [1] 7
# # >
# #     > nrow(table_s3_grubman_inh_neg )
# # > [1] 2
# #
# # > nrow(table_s3_grubman_oligo_neg)
# # > [1] 4
# #
# # > nrow(table_s3_grubman_OPC_neg)
# # [1] 1
# #---------------
# table_s3_grubman_astro_neg <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "astro" & grubman_logfc <= 0)
#
# table_s3_grubman_mg_neg <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "mg"& grubman_logfc <= 0)
#
# table_s3_grubman_ex_neg  <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "neuron (excitatory)" & grubman_logfc <= 0)
#
# table_s3_grubman_inh_neg <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "neuron (inhibitory)" & grubman_logfc <= 0)
#
# table_s3_grubman_oligo_neg <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "oligo" & grubman_logfc <= 0)
#
# table_s3_grubman_OPC_neg <- table_s3_grubman|>
#     dplyr::filter(`cell type` == "OPC" & grubman_logfc <= 0)
#
#
#
# df_mathys_list_directions <- list(
#     table_s3_grubman_astro_pos, table_s3_grubman_mg_pos, table_s3_grubman_ex_pos, table_s3_grubman_inh_pos,
#     table_s3_grubman_oligo_pos, table_s3_grubman_OPC_pos, table_s3_grubman_astro_neg,
#     table_s3_grubman_mg_neg, table_s3_grubman_ex_neg, table_s3_grubman_inh_neg,
#     table_s3_grubman_oligo_neg, table_s3_grubman_OPC_neg
# )
#
# res_2 <- purrr::map(df_mathys_list_directions, get_ensembl, Genes, "Genes")
#
#
# grubman_directions_geneList <- list(
#     grubman_mathys_astro_pos = res_2[[1]]$gene_ensembl_id,
#     grubman_mathys_mg_pos = res_2[[2]]$gene_ensembl_id,
#     grubman_mathys_ex_pos = res_2[[3]]$gene_ensembl_id,
#     grubman_mathys_inh_pos = res_2[[4]]$gene_ensembl_id,
#     grubman_mathys_oligo_pos = res_2[[5]]$gene_ensembl_id,
#     grubman_mathys_OPC_pos = res_2[[6]]$gene_ensembl_id,
#     grubman_mathys_astro_neg = res_2[[7]]$gene_ensembl_id,
#     grubman_mathys_mg_neg = res_2[[8]]$gene_ensembl_id,
#     grubman_mathys_ex_neg = res_2[[9]]$gene_ensembl_id,
#     grubman_mathys_inh_neg = res_2[[10]]$gene_ensembl_id,
#     grubman_mathys_oligo_neg = res_2[[11]]$gene_ensembl_id,
#     grubman_mathys_OPC_neg = res_2[[12]]$gene_ensembl_id
#
# )
#
#
# grubman_directions_enrichment <- gene_set_enrichment(
#     grubman_directions_geneList,
#     fdr_cut = 0.1,
#     modeling_results = modeling_results,
#     model_type = "enrichment"
# )
#
# grubman_directions_depleted <- gene_set_enrichment(
#     grubman_directions_geneList,
#     fdr_cut = 0.1,
#     modeling_results = modeling_results,
#     model_type = "enrichment",
#     reverse = TRUE
# )
#
#
# ##### enrichment plotting #####
# output_dir <- here("plots", "14_external_gene_sets")
# pdf(paste0(output_dir, "/04_grubman_directions_enriched.pdf"), width = 15)
# gene_set_enrichment_plot(
#     grubman_directions_enrichment,
#     xlabs = unique(grubman_directions_enrichment$ID),
#     PThresh = 12,
#     ORcut = 1.30103,
#     enrichOnly = FALSE,
#     layerHeights = c(0, seq_len(length(unique(grubman_directions_enrichment$test)))) * 15,
#     mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
#         9,
#         "YlOrRd"
#     )))(50)),
#     cex = 1.2
# )
#
# dev.off()
#
#
# pdf(paste0(output_dir, "/04_grubman_directions_depleted.pdf"), width = 15)
# gene_set_enrichment_plot(
#     grubman_directions_depleted,
#     xlabs = unique(grubman_directions_depleted$ID),
#     PThresh = 12,
#     ORcut = 1.30103,
#     enrichOnly = FALSE,
#     layerHeights = c(0, seq_len(length(unique(grubman_directions_depleted$test)))) * 15,
#     mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
#         9,
#         "YlOrRd"
#     )))(50)),
#     cex = 1.2
# )
#
# dev.off()
#
#
#
# #### GENES WITH CONCORDANCE == TRUE  ####
# table_s3_grubman_astro_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "astro" & Concordance == TRUE)
#
# table_s3_grubman_mg_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "mg"& Concordance == TRUE)
#
# table_s3_grubman_ex_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "neuron (excitatory)" & Concordance == TRUE)
#
# table_s3_grubman_inh_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "neuron (inhibitory)" & Concordance == TRUE)
#
# table_s3_grubman_oligo_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "oligo" & Concordance == TRUE)
#
# table_s3_grubman_OPC_conc <- table_s3_grubman |>
#     dplyr::filter(`cell type` == "OPC" & Concordance == TRUE)
#
#
# df_mathys_concordance_list <- list(
#     table_s3_grubman_astro_conc, table_s3_grubman_mg_conc, table_s3_grubman_ex_conc,
#     table_s3_grubman_inh_conc, table_s3_grubman_oligo_conc, table_s3_grubman_OPC_conc
# )
#
# res_3 <- purrr::map(df_mathys_concordance_list, get_ensembl, Genes, "Genes")
#
# grubman_conc_geneList <- list(
#     grubman_mathys_astro_conc = res_3[[1]]$gene_ensembl_id,
#     grubman_mathys_mg_conc = res_3[[2]]$gene_ensembl_id,
#     grubman_mathys_ex_conc = res_3[[3]]$gene_ensembl_id,
#     grubman_mathys_inh_conc = res_3[[4]]$gene_ensembl_id,
#     grubman_mathys_oligo_conc = res_3[[5]]$gene_ensembl_id,
#     grubman_mathys_OPC_conc = res_3[[6]]$gene_ensembl_id
# )
#
#
# grubman_concordance_enrichment <- gene_set_enrichment(
#     grubman_conc_geneList,
#     fdr_cut = 0.1,
#     modeling_results = modeling_results,
#     model_type = "enrichment"
# )
#
# grubman_concordance_depleted <- gene_set_enrichment(
#     grubman_conc_geneList,
#     fdr_cut = 0.1,
#     modeling_results = modeling_results,
#     model_type = "enrichment",
#     reverse = TRUE
# )
#
# ##### enrichment plotting #####
# pdf(paste0(output_dir, "/04_grubman_concordance_enriched.pdf"), width = 15)
# gene_set_enrichment_plot(
#     grubman_concordance_enrichment,
#     xlabs = unique(grubman_concordance_enrichment$ID),
#     PThresh = 12,
#     ORcut = 1.30103,
#     enrichOnly = FALSE,
#     layerHeights = c(0, seq_len(length(unique(grubman_concordance_enrichment$test)))) * 15,
#     mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
#         9,
#         "YlOrRd"
#     )))(50)),
#     cex = 1.2
# )
#
# dev.off()
#
#
# pdf(paste0(output_dir, "/04_grubman_concordance_depleted.pdf"), width = 15)
# gene_set_enrichment_plot(
#     grubman_concordance_depleted,
#     xlabs = unique(grubman_concordance_depleted$ID),
#     PThresh = 12,
#     ORcut = 1.30103,
#     enrichOnly = FALSE,
#     layerHeights = c(0, seq_len(length(unique(grubman_concordance_depleted$test)))) * 15,
#     mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
#         9,
#         "YlOrRd"
#     )))(50)),
#     cex = 1.2
# )
#
# dev.off()
#
#
#


# ####### Reproducibility information #####

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
