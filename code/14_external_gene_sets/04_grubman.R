#### load relevant packages ####

library("sgejobs")

# sgejobs::job_single(
#     "grubman",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 04_grubman.R"
# )

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
#  For each cell type, use all genes (ignore Concordance), and separate by direction.

input_data <- "raw-data/GeneSets/2_snRNA-seq/2_Grubman et al_Entorhinal/Grubman et al.xlsx"

# comparison with Mathys
table_s3_grub <- read_excel(input_data,
    sheet = "Supplementary Table 3", skip = 10,
    col_names = TRUE
)
# > nrow(table_s3_grub)
# [1] 94                                                     col_names = TRUE)

# > colnames(table_s3_grub)
# [1] "...1"                                    "Genes"
# [3] "Concordance"                             "Grubman.LogFC (AD vs Control a priori)"
# [5] "Mathys.LogFC(no pathology vs pathology)" "cell type"


# > unique(table_s3_grub$`cell type`)
# [1] "astro"               "mg"                  "neuron (excitatory)"
# [4] "neuron (inhibitory)" "oligo"               "OPC"


# "Legend:
# cond = a priori condition (AD or Control)
# cellType = identified cell type using BRETIGEA markers
# patient = patient for which gene is associated with
# geneID = gene symbol
# logFC = log2FoldChange
# FDR = false discovery rate using Benajmini-Hochberg Method   "


# #microglia
# table_s4a <-read_excel(input_data, sheet = "Supplementary Table 4a", skip = 6,
# col_names = TRUE)
# # nrow(table_s4a)
# # [1] 146
# nrow(table_s4a |> dplyr::filter(FDR < 0.1))
# # > colnames(table_s4a)
# # [1] "...1"     "cond"     "cellType" "patient"  "geneID"   "logFC"    "FDR"
#
# #astrocytes
# table_s4b <-read_excel(input_data, sheet = "Supplementary Table 4b", skip = 6,
#                        col_names = TRUE)
# # > nrow(table_s4b)
# # [1] 1503
# nrow(table_s4b |> dplyr::filter(FDR < 0.1))
#
# #neurons
# table_s4c <-read_excel(input_data, sheet = "Supplementary Table 4c", skip = 6,
#                        col_names = TRUE)
# nrow(table_s4c |> dplyr::filter(FDR < 0.1))
#
# # > nrow(table_s4c)
# # [1] 496
# # > unique(table_s4c$cellType)
# # [1] "neuron"
# #neurons aren't divided into excitatory and inhibitory
#
# #oligos
# table_s4d <-read_excel(input_data, sheet = "Supplementary Table 4d", skip = 6,
#                        col_names = TRUE)
# # > nrow(table_s4d)
# # [1] 2076
#
# nrow(table_s4d |> dplyr::filter(FDR < 0.1))
# #opc
# table_s4e <-read_excel(input_data, sheet = "Supplementary Table 4e", skip = 6,
#                        col_names = TRUE)
#
# nrow(table_s4e |> dplyr::filter(FDR < 0.1))
# # > nrow(table_s4e)
# # [1] 379

## Already filtered for FRD
# > nrow(table_s4a |> dplyr::filter(FDR < 0.1))
# [1] 146
# > nrow(table_s4b |> dplyr::filter(FDR < 0.1))
# [1] 1503
# > nrow(table_s4c |> dplyr::filter(FDR < 0.1))
# [1] 496
# > nrow(table_s4d |> dplyr::filter(FDR < 0.1))
# [1] 2076
# > nrow(table_s4e |> dplyr::filter(FDR < 0.1))
# [1] 379


# > unique(table_s3_grub$`cell type`)
# [1] "astro"               "mg"                  "neuron (excitatory)"
# [4] "neuron (inhibitory)" "oligo"               "OPC"
table_s3_grubman_astro <- table_s3_grub |> dplyr::filter(`cell type` == "astro")
table_s3_grubman_mg <- table_s3_grub |> dplyr::filter(`cell type` == "mg")
table_s3_grubman_ex <- table_s3_grub |> dplyr::filter(`cell type` == "neuron (excitatory)")
table_s3_grubman_inh <- table_s3_grub |> dplyr::filter(`cell type` == "neuron (inhibitory)")
table_s3_grubman_oligo <- table_s3_grub |> dplyr::filter(`cell type` == "oligo")
table_s3_grubman_OPC <- table_s3_grub |> dplyr::filter(`cell type` == "OPC")

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
