
# library("sgejobs")

# sgejobs::job_single(
#     "leung",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 05_leung.R"
# )

#### load relevant packages ####

library("readxl")
library("spatialLIBD")
library("dplyr")
library("sessioninfo")
library("here")
library("scran")
library("purrr")

here()
### load get_ensemble function

source(here("code/14_external_gene_sets/get_ensembl_function.R"))

### load modeling results
load(here(
    "processed-data", "11_grey_matter_only", "wholegenome",
    "Visium_IF_AD_modeling_results.Rdata"
))

# Number of sets: 2 + 9 subpopulations * direction = 2 + 18 = 20 sets
# Note: snRNA-seq
# Table S1: try it with and without direction
# Direction available: logFC
# Statistics available: FDR (filtered to 10%)
# Table S2: split by subpopulation, try it with and without direction
# Direction available: logFC
# Statistics available: globalFDR (filtered to 10%)
# Ignore the comparison for now since some of them are very small


# braak comparison's of interest
# 197:371
# 474:552


input_dir <- here("raw-data", "GeneSets", "2_snRNA-seq",
                  "3_Leung et al", "Leung et al.xlsx")

leung_1 <- read_excel(input_dir,
    sheet = "Supplementary Table 1", col_names = TRUE
)
# nrow(leung_1)
# # [1] 434
leung_1 <- leung_1 |> dplyr::filter(FDR < 0.1)
# nrow(leung_1)
# # [1] 434

leung_1_up <- leung_1 |> dplyr::filter(logFC > 0)
# nrow(leung_1_up)
# [1] 195

leung_1_down <- leung_1 |> dplyr::filter(logFC <= 0)
# nrow(leung_1_down)
# 239


leung_2 <- read_excel(input_dir,
    sheet = "Supplementary Table 2", col_names = TRUE
)
leung_2 <- leung_2 |> dplyr::filter(globalFDR < 0.1)
nrow(leung_2)
# 1117
unique(leung_2$subpopulation)
# > unique(leung_2$subpopulation)
# [1] "EC:Exc.s0" "EC:Exc.s1" "EC:Exc.s2" "EC:Exc.s3" "EC:Exc.s4"
# [6] "EC:Exc.s5" "EC:Exc.s6" "EC:Exc.s7" "EC:Exc.s8"

leung_2_up <- leung_2 |> dplyr::filter(logFC > 0)
nrow(leung_2_up)
# > nrow(leung_2_up)
# [1] 614
leung_2_down <- leung_2 |> dplyr::filter(logFC <= 0)
nrow(leung_2_down)
# [1] 503


#### Sort gene sets ####

# leung_2_s0 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s0")
# leung_2_s1 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s1")
# leung_2_s2 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s2")
# leung_2_s3 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s3")
# leung_2_s4 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s4")
# leung_2_s5 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s5")
# leung_2_s6 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s6")
# leung_2_s7 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s7")
# leung_2_s8 <- leung_2 |> dplyr::filter(subpopulation == "EC:Exc.s8")


leung_2_up_s0 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s0") # 58
leung_2_up_s1 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s1") # 33
leung_2_up_s2 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s2") # 79
leung_2_up_s3 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s3") # 60
leung_2_up_s4 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s4") # 48
leung_2_up_s5 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s5") # 108
leung_2_up_s6 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s6") # 48
leung_2_up_s7 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s7") # 73
leung_2_up_s8 <- leung_2_up |> dplyr::filter(subpopulation == "EC:Exc.s8") # 107


leung_2_down_s0 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s0") # 51
leung_2_down_s1 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s1") # 49
leung_2_down_s2 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s2") # 108
leung_2_down_s3 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s3") # 24
leung_2_down_s4 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s4") # 41
leung_2_down_s5 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s5") # 64
leung_2_down_s6 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s6") # 35
leung_2_down_s7 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s7") # 61
leung_2_down_s8 <- leung_2_down |> dplyr::filter(subpopulation == "EC:Exc.s8") # 70

#### Create gene ensembl sets ####
df_list <- list(
    leung_1_up, leung_1_down, leung_2_up_s0, leung_2_up_s1, leung_2_up_s2,
    leung_2_up_s3, leung_2_up_s4, leung_2_up_s5, leung_2_up_s6, leung_2_up_s7,
    leung_2_up_s8, leung_2_down_s0, leung_2_down_s1, leung_2_down_s2, leung_2_down_s3,
    leung_2_down_s4, leung_2_down_s5, leung_2_down_s6, leung_2_down_s7, leung_2_down_s8
)

res_1 <- purrr::map(df_list, get_ensembl, gene, "gene")



leung_geneList <- list(
    leung_1_up = res_1[[1]]$gene_ensembl_id,
    leung_1_down = res_1[[2]]$gene_ensembl_id,
    leung_2_up_s0 = res_1[[3]]$gene_ensembl_id,
    leung_2_up_s1 = res_1[[4]]$gene_ensembl_id,
    leung_2_up_s2 = res_1[[5]]$gene_ensembl_id,
    leung_2_up_s3 = res_1[[6]]$gene_ensembl_id,
    leung_2_up_s4 = res_1[[7]]$gene_ensembl_id,
    leung_2_up_s5 = res_1[[8]]$gene_ensembl_id,
    leung_2_up_s6 = res_1[[9]]$gene_ensembl_id,
    leung_2_up_s7 = res_1[[10]]$gene_ensembl_id,
    leung_2_up_s8 = res_1[[11]]$gene_ensembl_id,
    leung_2_down_s0 = res_1[[12]]$gene_ensembl_id,
    leung_2_down_s1 = res_1[[13]]$gene_ensembl_id,
    leung_2_down_s2 = res_1[[14]]$gene_ensembl_id,
    leung_2_down_s3 = res_1[[15]]$gene_ensembl_id,
    leung_2_down_s4 = res_1[[16]]$gene_ensembl_id,
    leung_2_down_s5 = res_1[[17]]$gene_ensembl_id,
    leung_2_down_s6 = res_1[[18]]$gene_ensembl_id,
    leung_2_down_s7 = res_1[[19]]$gene_ensembl_id,
    leung_2_down_s8 = res_1[[20]]$gene_ensembl_id
)

#### Perform enrichment analysis ####
leung_enrichment <- gene_set_enrichment(
    leung_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

leung_enrichment
# OR        Pval      test              ID model_type fdr_cut
# 1    0.000000 1.000000000      none      leung_1_up enrichment     0.1
# 2    0.000000 1.000000000      none    leung_1_down enrichment     0.1
# 3    0.000000 1.000000000      none   leung_2_up_s0 enrichment     0.1
# 4    0.000000 1.000000000      none   leung_2_up_s1 enrichment     0.1
# 5    0.000000 1.000000000      none   leung_2_up_s2 enrichment     0.1
# 6    0.000000 1.000000000      none   leung_2_up_s3 enrichment     0.1
# 7    0.000000 1.000000000      none   leung_2_up_s4 enrichment     0.1
# 8    0.000000 1.000000000      none   leung_2_up_s5 enrichment     0.1
# 9    0.000000 1.000000000      none   leung_2_up_s6 enrichment     0.1
# 10   0.000000 1.000000000      none   leung_2_up_s7 enrichment     0.1
# 11   0.000000 1.000000000      none   leung_2_up_s8 enrichment     0.1
# 12   0.000000 1.000000000      none leung_2_down_s0 enrichment     0.1
# 13   0.000000 1.000000000      none leung_2_down_s1 enrichment     0.1
# 14   0.000000 1.000000000      none leung_2_down_s2 enrichment     0.1
# 15   0.000000 1.000000000      none leung_2_down_s3 enrichment     0.1
# 16   0.000000 1.000000000      none leung_2_down_s4 enrichment     0.1
# 17   0.000000 1.000000000      none leung_2_down_s5 enrichment     0.1
# 18   0.000000 1.000000000      none leung_2_down_s6 enrichment     0.1
# 19   0.000000 1.000000000      none leung_2_down_s7 enrichment     0.1
# 20   0.000000 1.000000000      none leung_2_down_s8 enrichment     0.1
# 21   0.000000 1.000000000       Ab+      leung_1_up enrichment     0.1
# 22   0.000000 1.000000000       Ab+    leung_1_down enrichment     0.1
# 23   0.000000 1.000000000       Ab+   leung_2_up_s0 enrichment     0.1
# 24   0.000000 1.000000000       Ab+   leung_2_up_s1 enrichment     0.1
# 25   0.000000 1.000000000       Ab+   leung_2_up_s2 enrichment     0.1
# 26   0.000000 1.000000000       Ab+   leung_2_up_s3 enrichment     0.1
# 27   0.000000 1.000000000       Ab+   leung_2_up_s4 enrichment     0.1
# 28   0.000000 1.000000000       Ab+   leung_2_up_s5 enrichment     0.1
# 29   0.000000 1.000000000       Ab+   leung_2_up_s6 enrichment     0.1
# 30   0.000000 1.000000000       Ab+   leung_2_up_s7 enrichment     0.1
# 31   0.000000 1.000000000       Ab+   leung_2_up_s8 enrichment     0.1
# 32   0.000000 1.000000000       Ab+ leung_2_down_s0 enrichment     0.1
# 33   0.000000 1.000000000       Ab+ leung_2_down_s1 enrichment     0.1
# 34   0.000000 1.000000000       Ab+ leung_2_down_s2 enrichment     0.1
# 35   0.000000 1.000000000       Ab+ leung_2_down_s3 enrichment     0.1
# 36   0.000000 1.000000000       Ab+ leung_2_down_s4 enrichment     0.1
# 37   0.000000 1.000000000       Ab+ leung_2_down_s5 enrichment     0.1
# 38   0.000000 1.000000000       Ab+ leung_2_down_s6 enrichment     0.1
# 39   0.000000 1.000000000       Ab+ leung_2_down_s7 enrichment     0.1
# 40   0.000000 1.000000000       Ab+ leung_2_down_s8 enrichment     0.1
# 41   0.000000 1.000000000  next_Ab+      leung_1_up enrichment     0.1
# 42   1.161948 0.582638915  next_Ab+    leung_1_down enrichment     0.1
# 43   0.000000 1.000000000  next_Ab+   leung_2_up_s0 enrichment     0.1
# 44   0.000000 1.000000000  next_Ab+   leung_2_up_s1 enrichment     0.1
# 45   0.000000 1.000000000  next_Ab+   leung_2_up_s2 enrichment     0.1
# 46   0.000000 1.000000000  next_Ab+   leung_2_up_s3 enrichment     0.1
# 47   0.000000 1.000000000  next_Ab+   leung_2_up_s4 enrichment     0.1
# 48  14.902987 0.001509311  next_Ab+   leung_2_up_s5 enrichment     0.1
# 49   0.000000 1.000000000  next_Ab+   leung_2_up_s6 enrichment     0.1
# 50   8.365086 0.119606858  next_Ab+   leung_2_up_s7 enrichment     0.1
# 51  13.951439 0.075794583  next_Ab+   leung_2_up_s8 enrichment     0.1
# 52   0.000000 1.000000000  next_Ab+ leung_2_down_s0 enrichment     0.1
# 53   0.000000 1.000000000  next_Ab+ leung_2_down_s1 enrichment     0.1
# 54   0.000000 1.000000000  next_Ab+ leung_2_down_s2 enrichment     0.1
# 55   0.000000 1.000000000  next_Ab+ leung_2_down_s3 enrichment     0.1
# 56   0.000000 1.000000000  next_Ab+ leung_2_down_s4 enrichment     0.1
# 57   0.000000 1.000000000  next_Ab+ leung_2_down_s5 enrichment     0.1
# 58   0.000000 1.000000000  next_Ab+ leung_2_down_s6 enrichment     0.1
# 59   0.000000 1.000000000  next_Ab+ leung_2_down_s7 enrichment     0.1
# 60   0.000000 1.000000000  next_Ab+ leung_2_down_s8 enrichment     0.1
# 61   0.000000 1.000000000       pT+      leung_1_up enrichment     0.1
# 62   0.000000 1.000000000       pT+    leung_1_down enrichment     0.1
# 63   0.000000 1.000000000       pT+   leung_2_up_s0 enrichment     0.1
# 64   0.000000 1.000000000       pT+   leung_2_up_s1 enrichment     0.1
# 65   0.000000 1.000000000       pT+   leung_2_up_s2 enrichment     0.1
# 66   0.000000 1.000000000       pT+   leung_2_up_s3 enrichment     0.1
# 67   0.000000 1.000000000       pT+   leung_2_up_s4 enrichment     0.1
# 68   0.000000 1.000000000       pT+   leung_2_up_s5 enrichment     0.1
# 69   0.000000 1.000000000       pT+   leung_2_up_s6 enrichment     0.1
# 70   0.000000 1.000000000       pT+   leung_2_up_s7 enrichment     0.1
# 71   0.000000 1.000000000       pT+   leung_2_up_s8 enrichment     0.1
# 72   0.000000 1.000000000       pT+ leung_2_down_s0 enrichment     0.1
# 73   0.000000 1.000000000       pT+ leung_2_down_s1 enrichment     0.1
# 74   0.000000 1.000000000       pT+ leung_2_down_s2 enrichment     0.1
# 75   0.000000 1.000000000       pT+ leung_2_down_s3 enrichment     0.1
# 76   0.000000 1.000000000       pT+ leung_2_down_s4 enrichment     0.1
# 77   0.000000 1.000000000       pT+ leung_2_down_s5 enrichment     0.1
# 78   0.000000 1.000000000       pT+ leung_2_down_s6 enrichment     0.1
# 79   0.000000 1.000000000       pT+ leung_2_down_s7 enrichment     0.1
# 80   0.000000 1.000000000       pT+ leung_2_down_s8 enrichment     0.1
# 81   0.000000 1.000000000  next_pT+      leung_1_up enrichment     0.1
# 82   0.000000 1.000000000  next_pT+    leung_1_down enrichment     0.1
# 83   0.000000 1.000000000  next_pT+   leung_2_up_s0 enrichment     0.1
# 84   0.000000 1.000000000  next_pT+   leung_2_up_s1 enrichment     0.1
# 85   0.000000 1.000000000  next_pT+   leung_2_up_s2 enrichment     0.1
# 86   0.000000 1.000000000  next_pT+   leung_2_up_s3 enrichment     0.1
# 87   0.000000 1.000000000  next_pT+   leung_2_up_s4 enrichment     0.1
# 88   0.000000 1.000000000  next_pT+   leung_2_up_s5 enrichment     0.1
# 89   0.000000 1.000000000  next_pT+   leung_2_up_s6 enrichment     0.1
# 90   0.000000 1.000000000  next_pT+   leung_2_up_s7 enrichment     0.1
# 91   0.000000 1.000000000  next_pT+   leung_2_up_s8 enrichment     0.1
# 92   0.000000 1.000000000  next_pT+ leung_2_down_s0 enrichment     0.1
# 93   0.000000 1.000000000  next_pT+ leung_2_down_s1 enrichment     0.1
# 94   0.000000 1.000000000  next_pT+ leung_2_down_s2 enrichment     0.1
# 95   0.000000 1.000000000  next_pT+ leung_2_down_s3 enrichment     0.1
# 96   0.000000 1.000000000  next_pT+ leung_2_down_s4 enrichment     0.1
# 97   0.000000 1.000000000  next_pT+ leung_2_down_s5 enrichment     0.1
# 98   0.000000 1.000000000  next_pT+ leung_2_down_s6 enrichment     0.1
# 99   0.000000 1.000000000  next_pT+ leung_2_down_s7 enrichment     0.1
# 100  0.000000 1.000000000  next_pT+ leung_2_down_s8 enrichment     0.1
# 101  0.000000 1.000000000      both      leung_1_up enrichment     0.1
# 102  0.000000 1.000000000      both    leung_1_down enrichment     0.1
# 103  0.000000 1.000000000      both   leung_2_up_s0 enrichment     0.1
# 104  0.000000 1.000000000      both   leung_2_up_s1 enrichment     0.1
# 105  0.000000 1.000000000      both   leung_2_up_s2 enrichment     0.1
# 106  0.000000 1.000000000      both   leung_2_up_s3 enrichment     0.1
# 107  0.000000 1.000000000      both   leung_2_up_s4 enrichment     0.1
# 108  0.000000 1.000000000      both   leung_2_up_s5 enrichment     0.1
# 109  0.000000 1.000000000      both   leung_2_up_s6 enrichment     0.1
# 110  0.000000 1.000000000      both   leung_2_up_s7 enrichment     0.1
# 111  0.000000 1.000000000      both   leung_2_up_s8 enrichment     0.1
# 112  0.000000 1.000000000      both leung_2_down_s0 enrichment     0.1
# 113  0.000000 1.000000000      both leung_2_down_s1 enrichment     0.1
# 114  0.000000 1.000000000      both leung_2_down_s2 enrichment     0.1
# 115  0.000000 1.000000000      both leung_2_down_s3 enrichment     0.1
# 116  0.000000 1.000000000      both leung_2_down_s4 enrichment     0.1
# 117  0.000000 1.000000000      both leung_2_down_s5 enrichment     0.1
# 118  0.000000 1.000000000      both leung_2_down_s6 enrichment     0.1
# 119  0.000000 1.000000000      both leung_2_down_s7 enrichment     0.1
# 120  0.000000 1.000000000      both leung_2_down_s8 enrichment     0.1
# 121  0.000000 1.000000000 next_both      leung_1_up enrichment     0.1
# 122  0.000000 1.000000000 next_both    leung_1_down enrichment     0.1
# 123  0.000000 1.000000000 next_both   leung_2_up_s0 enrichment     0.1
# 124  0.000000 1.000000000 next_both   leung_2_up_s1 enrichment     0.1
# 125  0.000000 1.000000000 next_both   leung_2_up_s2 enrichment     0.1
# 126  0.000000 1.000000000 next_both   leung_2_up_s3 enrichment     0.1
# 127  0.000000 1.000000000 next_both   leung_2_up_s4 enrichment     0.1
# 128  0.000000 1.000000000 next_both   leung_2_up_s5 enrichment     0.1
# 129  0.000000 1.000000000 next_both   leung_2_up_s6 enrichment     0.1
# 130  0.000000 1.000000000 next_both   leung_2_up_s7 enrichment     0.1
# 131  0.000000 1.000000000 next_both   leung_2_up_s8 enrichment     0.1
# 132  0.000000 1.000000000 next_both leung_2_down_s0 enrichment     0.1
# 133  0.000000 1.000000000 next_both leung_2_down_s1 enrichment     0.1
# 134  0.000000 1.000000000 next_both leung_2_down_s2 enrichment     0.1
# 135  0.000000 1.000000000 next_both leung_2_down_s3 enrichment     0.1
# 136  0.000000 1.000000000 next_both leung_2_down_s4 enrichment     0.1
# 137  0.000000 1.000000000 next_both leung_2_down_s5 enrichment     0.1
# 138  0.000000 1.000000000 next_both leung_2_down_s6 enrichment     0.1
# 139  0.000000 1.000000000 next_both leung_2_down_s7 enrichment     0.1
# 140  0.000000 1.000000000 next_both leung_2_down_s8 enrichment     0.1



leung_depleted <- gene_set_enrichment(
    leung_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE
)


##### enrichment plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/05_leung_enriched.pdf"), width = 15)
gene_set_enrichment_plot(
    leung_enrichment,
    xlabs = unique(leung_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(leung_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        9,
        "YlOrRd"
    )))(50)),
    cex = 1.2
)

dev.off()

pdf(paste0(output_dir, "/05_leung_depleted.pdf"), width = 15)
gene_set_enrichment_plot(
    leung_depleted,
    xlabs = unique(leung_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(leung_depleted$test)))) * 15,
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
