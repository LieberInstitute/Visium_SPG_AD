

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
table_s3_grubman$grubman_logfc <- table_s3_grub$`Grubman.LogFC (AD vs Control a priori)`
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
    grubman_mathys_astro = res_1[[1]]$gene_ensembl_id,
    grubman_mathys_mg = res_1[[2]]$gene_ensembl_id,
    grubman_mathys_ex = res_1[[3]]$gene_ensembl_id,
    grubman_mathys_inh = res_1[[4]]$gene_ensembl_id,
    grubman_mathys_oligo = res_1[[5]]$gene_ensembl_id,
    grubman_mathys_OPC = res_1[[6]]$gene_ensembl_id
)


grubman_enrichment <- gene_set_enrichment(
    grubman_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

grubman_depleted <- gene_set_enrichment(
    grubman_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)

#### grubman_enrichment ####
# > grubman_enrichment
# OR       Pval      test                  ID model_type fdr_cut
# 1   0.000000 1.00000000      none grub_v_mathys_astro enrichment     0.1
# 2   0.000000 1.00000000      none    grub_v_mathys_mg enrichment     0.1
# 3   0.000000 1.00000000      none    grub_v_mathys_ex enrichment     0.1
# 4   0.000000 1.00000000      none   grub_v_mathys_inh enrichment     0.1
# 5   0.000000 1.00000000      none grub_v_mathys_oligo enrichment     0.1
# 6   0.000000 1.00000000      none   grub_v_mathys_OPC enrichment     0.1
# 7   0.000000 1.00000000      none       table_4_astro enrichment     0.1
# 8   0.000000 1.00000000      none          table_4_mg enrichment     0.1
# 9   0.000000 1.00000000      none     table_4_neurons enrichment     0.1
# 10  0.000000 1.00000000      none         table_4_oli enrichment     0.1
# 11  0.000000 1.00000000      none         table_4_OPC enrichment     0.1
# 12  0.000000 1.00000000       Ab+ grub_v_mathys_astro enrichment     0.1
# 13  0.000000 1.00000000       Ab+    grub_v_mathys_mg enrichment     0.1
# 14  0.000000 1.00000000       Ab+    grub_v_mathys_ex enrichment     0.1
# 15  0.000000 1.00000000       Ab+   grub_v_mathys_inh enrichment     0.1
# 16  0.000000 1.00000000       Ab+ grub_v_mathys_oligo enrichment     0.1
# 17  0.000000 1.00000000       Ab+   grub_v_mathys_OPC enrichment     0.1
# 18  1.427570 0.53050297       Ab+       table_4_astro enrichment     0.1
# 19  0.000000 1.00000000       Ab+          table_4_mg enrichment     0.1
# 20  0.000000 1.00000000       Ab+     table_4_neurons enrichment     0.1
# 21  1.363989 0.54606201       Ab+         table_4_oli enrichment     0.1
# 22  0.000000 1.00000000       Ab+         table_4_OPC enrichment     0.1
# 23  6.689132 0.14594855  next_Ab+ grub_v_mathys_astro enrichment     0.1
# 24 23.886371 0.04733386  next_Ab+    grub_v_mathys_mg enrichment     0.1
# 25  9.297019 0.10884859  next_Ab+    grub_v_mathys_ex enrichment     0.1
# 26  0.000000 1.00000000  next_Ab+   grub_v_mathys_inh enrichment     0.1
# 27  0.000000 1.00000000  next_Ab+ grub_v_mathys_oligo enrichment     0.1
# 28  0.000000 1.00000000  next_Ab+   grub_v_mathys_OPC enrichment     0.1
# 29  1.242637 0.60269755  next_Ab+       table_4_astro enrichment     0.1
# 30  2.159499 0.37796767  next_Ab+          table_4_mg enrichment     0.1
# 31  2.018614 0.19912936  next_Ab+     table_4_neurons enrichment     0.1
# 32  1.458234 0.44076336  next_Ab+         table_4_oli enrichment     0.1
# 33  2.668662 0.11313297  next_Ab+         table_4_OPC enrichment     0.1
# 34  0.000000 1.00000000       pT+ grub_v_mathys_astro enrichment     0.1
# 35  0.000000 1.00000000       pT+    grub_v_mathys_mg enrichment     0.1
# 36  0.000000 1.00000000       pT+    grub_v_mathys_ex enrichment     0.1
# 37  0.000000 1.00000000       pT+   grub_v_mathys_inh enrichment     0.1
# 38  0.000000 1.00000000       pT+ grub_v_mathys_oligo enrichment     0.1
# 39  0.000000 1.00000000       pT+   grub_v_mathys_OPC enrichment     0.1
# 40  0.000000 1.00000000       pT+       table_4_astro enrichment     0.1
# 41  0.000000 1.00000000       pT+          table_4_mg enrichment     0.1
# 42  0.000000 1.00000000       pT+     table_4_neurons enrichment     0.1
# 43  0.000000 1.00000000       pT+         table_4_oli enrichment     0.1
# 44  0.000000 1.00000000       pT+         table_4_OPC enrichment     0.1
# 45  0.000000 1.00000000  next_pT+ grub_v_mathys_astro enrichment     0.1
# 46  0.000000 1.00000000  next_pT+    grub_v_mathys_mg enrichment     0.1
# 47  0.000000 1.00000000  next_pT+    grub_v_mathys_ex enrichment     0.1
# 48  0.000000 1.00000000  next_pT+   grub_v_mathys_inh enrichment     0.1
# 49  0.000000 1.00000000  next_pT+ grub_v_mathys_oligo enrichment     0.1
# 50  0.000000 1.00000000  next_pT+   grub_v_mathys_OPC enrichment     0.1
# 51  0.000000 1.00000000  next_pT+       table_4_astro enrichment     0.1
# 52  0.000000 1.00000000  next_pT+          table_4_mg enrichment     0.1
# 53  0.000000 1.00000000  next_pT+     table_4_neurons enrichment     0.1
# 54  0.000000 1.00000000  next_pT+         table_4_oli enrichment     0.1
# 55  0.000000 1.00000000  next_pT+         table_4_OPC enrichment     0.1
# 56  0.000000 1.00000000      both grub_v_mathys_astro enrichment     0.1
# 57  0.000000 1.00000000      both    grub_v_mathys_mg enrichment     0.1
# 58  0.000000 1.00000000      both    grub_v_mathys_ex enrichment     0.1
# 59  0.000000 1.00000000      both   grub_v_mathys_inh enrichment     0.1
# 60  0.000000 1.00000000      both grub_v_mathys_oligo enrichment     0.1
# 61  0.000000 1.00000000      both   grub_v_mathys_OPC enrichment     0.1
# 62  0.000000 1.00000000      both       table_4_astro enrichment     0.1
# 63  0.000000 1.00000000      both          table_4_mg enrichment     0.1
# 64  0.000000 1.00000000      both     table_4_neurons enrichment     0.1
# 65  0.000000 1.00000000      both         table_4_oli enrichment     0.1
# 66  0.000000 1.00000000      both         table_4_OPC enrichment     0.1
# 67  0.000000 1.00000000 next_both grub_v_mathys_astro enrichment     0.1
# 68  0.000000 1.00000000 next_both    grub_v_mathys_mg enrichment     0.1
# 69  0.000000 1.00000000 next_both    grub_v_mathys_ex enrichment     0.1
# 70  0.000000 1.00000000 next_both   grub_v_mathys_inh enrichment     0.1
# 71  0.000000 1.00000000 next_both grub_v_mathys_oligo enrichment     0.1
# 72  0.000000 1.00000000 next_both   grub_v_mathys_OPC enrichment     0.1
# 73  0.000000 1.00000000 next_both       table_4_astro enrichment     0.1
# 74  0.000000 1.00000000 next_both          table_4_mg enrichment     0.1
# 75  0.000000 1.00000000 next_both     table_4_neurons enrichment     0.1
# 76  0.000000 1.00000000 next_both         table_4_oli enrichment     0.1
# 77  0.000000 1.00000000 next_both         table_4_OPC enrichment     0.1
#### ####

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


#### GENES SPLIT BY DIRECTION  ####

table_s3_grubman_astro_pos <- table_s3_grubman|>
    dplyr::filter(`cell type` == "astro" & grubman_logfc > 0)

table_s3_grubman_mg_pos <- table_s3_grubman|>
    dplyr::filter(`cell type` == "mg"& grubman_logfc > 0)

table_s3_grubman_ex_pos  <- table_s3_grubman|>
    dplyr::filter(`cell type` == "neuron (excitatory)" & grubman_logfc > 0)

table_s3_grubman_inh_pos <- table_s3_grubman|>
    dplyr::filter(`cell type` == "neuron (inhibitory)" & grubman_logfc > 0)

table_s3_grubman_oligo_pos <- table_s3_grubman|>
    dplyr::filter(`cell type` == "oligo" & grubman_logfc > 0)

table_s3_grubman_OPC_pos <- table_s3_grubman|>
    dplyr::filter(`cell type` == "OPC" & grubman_logfc > 0)

#---------------
table_s3_grubman_astro_neg <- table_s3_grubman|>
    dplyr::filter(`cell type` == "astro" & grubman_logfc <= 0)

table_s3_grubman_mg_neg <- table_s3_grubman|>
    dplyr::filter(`cell type` == "mg"& grubman_logfc <= 0)

table_s3_grubman_ex_neg  <- table_s3_grubman|>
    dplyr::filter(`cell type` == "neuron (excitatory)" & grubman_logfc <= 0)

table_s3_grubman_inh_neg <- table_s3_grubman|>
    dplyr::filter(`cell type` == "neuron (inhibitory)" & grubman_logfc <= 0)

table_s3_grubman_oligo_neg <- table_s3_grubman|>
    dplyr::filter(`cell type` == "oligo" & grubman_logfc <= 0)

table_s3_grubman_OPC_neg <- table_s3_grubman|>
    dplyr::filter(`cell type` == "OPC" & grubman_logfc <= 0)



df_mathys_list_directions <- list(
    table_s3_grubman_astro_pos, table_s3_grubman_mg_pos, table_s3_grubman_ex_pos, table_s3_grubman_inh_pos,
    table_s3_grubman_oligo_pos, table_s3_grubman_OPC_pos, table_s3_grubman_astro_neg,
    table_s3_grubman_mg_neg, table_s3_grubman_ex_neg, table_s3_grubman_inh_neg,
    table_s3_grubman_oligo_neg, table_s3_grubman_OPC_neg
)

res_2 <- purrr::map(df_mathys_list_directions, get_ensembl, Genes, "Genes")


grubman_directions_geneList <- list(
    grubman_mathys_astro_pos = res_2[[1]]$gene_ensembl_id,
    grubman_mathys_mg_pos = res_2[[2]]$gene_ensembl_id,
    grubman_mathys_ex_pos = res_2[[3]]$gene_ensembl_id,
    grubman_mathys_inh_pos = res_2[[4]]$gene_ensembl_id,
    grubman_mathys_oligo_pos = res_2[[5]]$gene_ensembl_id,
    grubman_mathys_OPC_pos = res_2[[6]]$gene_ensembl_id,
    grubman_mathys_astro_neg = res_2[[7]]$gene_ensembl_id,
    grubman_mathys_mg_neg = res_2[[8]]$gene_ensembl_id,
    grubman_mathys_ex_neg = res_2[[9]]$gene_ensembl_id,
    grubman_mathys_inh_neg = res_2[[10]]$gene_ensembl_id,
    grubman_mathys_oligo_neg = res_2[[11]]$gene_ensembl_id,
    grubman_mathys_OPC_neg = res_2[[12]]$gene_ensembl_id

)


grubman_directions_enrichment <- gene_set_enrichment(
    grubman_directions_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

grubman_directions_depleted <- gene_set_enrichment(
    grubman_directions_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)


##### enrichment plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/04_grubman_directions_enriched.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_directions_enrichment,
    xlabs = unique(grubman_directions_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_directions_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()


pdf(paste0(output_dir, "/04_grubman_directions_depleted.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_directions_depleted,
    xlabs = unique(grubman_directions_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_directions_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()



#### GENES WITH CONCORDANCE == TRUE  ####
table_s3_grubman_astro_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "astro" & Concordance == TRUE)

table_s3_grubman_mg_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "mg"& Concordance == TRUE)

table_s3_grubman_ex_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "neuron (excitatory)" & Concordance == TRUE)

table_s3_grubman_inh_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "neuron (inhibitory)" & Concordance == TRUE)

table_s3_grubman_oligo_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "oligo" & Concordance == TRUE)

table_s3_grubman_OPC_conc <- table_s3_grubman |>
    dplyr::filter(`cell type` == "OPC" & Concordance == TRUE)


# > nrow(table_s3_grubman_astro_conc )
# [1] 32
# > nrow(table_s3_grubman_mg_conc)
# [1] 10
# > nrow(table_s3_grubman_ex_conc)
# [1] 20
# > nrow(table_s3_grubman_inh_conc)
# [1] 3
# > nrow(table_s3_grubman_oligo_conc)
# [1] 20
# > nrow(table_s3_grubman_OPC_conc)
# [1] 4


df_mathys_concordance_list <- list(
    table_s3_grubman_astro_conc, table_s3_grubman_mg_conc, table_s3_grubman_ex_conc,
    table_s3_grubman_inh_conc, table_s3_grubman_oligo_conc, table_s3_grubman_OPC_conc
)

res_3 <- purrr::map(df_mathys_concordance_list, get_ensembl, Genes, "Genes")

grubman_geneList <- list(
    grubman_mathys_astro_conc = res_3[[1]]$gene_ensembl_id,
    grubman_mathys_mg_conc = res_3[[2]]$gene_ensembl_id,
    grubman_mathys_ex_conc = res_3[[3]]$gene_ensembl_id,
    grubman_mathys_inh_conc = res_3[[4]]$gene_ensembl_id,
    grubman_mathys_oligo_conc = res_3[[5]]$gene_ensembl_id,
    grubman_mathys_OPC_conc = res_3[[6]]$gene_ensembl_id
)


grubman_concordance_enrichment <- gene_set_enrichment(
    grubman_conc_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

grubman_concordance_depleted <- gene_set_enrichment(
    grubman_conc_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)

##### enrichment plotting #####
pdf(paste0(output_dir, "/04_grubman_concordance_enriched.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_concordance_enrichment,
    xlabs = unique(grubman_concordance_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_concordance_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()


pdf(paste0(output_dir, "/04_grubman_concordance_depleted.pdf"), width = 15)
gene_set_enrichment_plot(
    grubman_concordance_depleted,
    xlabs = unique(grubman_concordance_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(grubman_concordance_depleted$test)))) * 15,
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
