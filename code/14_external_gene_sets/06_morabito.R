#load libraries

library('readxl')
library('spatialLIBD')
library('dplyr')
library('sessioninfo')
library('here')
library('scran')
library('purrr')


# Number of sets: 6 cell types * direction (3) * 2 assays (ATAC and RNA) = 36 sets
# Note: snATAC-seq and snRNA-seq
# (low priority since it's ATAC-seq with no assigned genes)

#  Table S1: split by celltype, try it with and without direction
# Direction available: avg_logFC
# Statistics available: p_val_adj (looks filtered already)

#Table S5: try it with and without direction
# Direction available: avg_logFC
# Statistics available: p_val_adj
# Genes are not assigned. We have the Peak coordinates. We could use GenomicRanges::findOverlaps() to compare against the annotation.
#  Let's try subsetting to Peaks that overlap a single gene first (including the UTRs). If they don't overlap genes or overlap to 2 or more genes, we'll ignore them. But it could be that we lose too many peaks that way.
#     Another option is to find the distance to the closest 5' of a gene. Maybe with a maximum distance filter.



### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensembl_function.R'))

### load modeling results
load(here('processed-data','11_grey_matter_only','wholegenome',
          'Visium_IF_AD_modeling_results.Rdata'))

table_1 <- read_excel("raw-data/GeneSets/3_snATAC-seq/Table S1_snRNAseq.xlsx",
                      sheet = "Supplementary Data 1e", col_names = TRUE, skip = 2)
head(table_1)
# > nrow(table_1)
# [1] 3085


# p_val avg_logFC pct.1 pct.2 p_val_adj gene    diff celltype
# <dbl>     <dbl> <dbl> <dbl>     <dbl> <chr>  <dbl> <chr>
# 1     0     0.717 0.603 0.231         0 XIST  0.372  ODC
# 2     0     0.652 0.448 0.254         0 PTPRM 0.194  ODC
# 3     0     0.492 0.96  0.922         0 NEAT1 0.0380 ODC
# 4     0     0.348 0.726 0.539         0 FKBP5 0.187  ODC
# 5     0     0.333 0.61  0.432         0 XPO1  0.178  ODC
# 6     0     0.328 0.379 0.212         0 CHOR… 0.167  ODC


#filter by FDR
table_1  <- table_1 |> dplyr::filter(p_val_adj < 0.1)



table_1_up <- table_1 |> dplyr::filter(avg_logFC > 0)
nrow(table_1_up)
# > nrow(table_1_up)
# [1]1495
table_1_down <- table_1 |> dplyr::filter(avg_logFC <= 0)
nrow(table_1_down)
# [1]  1590



#### Sort gene sets ####
unique(table_1$celltype)
#"ODC"     "MG"      "OPC"     "INH"     "EX"      "ASC"     "PER.END"

table_1_up_ODC <- table_1_up  |> dplyr::filter(celltype == "ODC") #327
table_1_up_MG <- table_1_up  |> dplyr::filter(celltype == "MG") #246
table_1_up_OPC <- table_1_up  |> dplyr::filter(celltype == "OPC") #175
table_1_up_INH <- table_1_up  |> dplyr::filter(celltype == "INH") #125
table_1_up_EX <- table_1_up  |> dplyr::filter(celltype == "EX") #231
table_1_up_ASC <- table_1_up  |> dplyr::filter(celltype == "ASC") #384
table_1_up_PER.END <- table_1_up  |> dplyr::filter(celltype == "PER.END") #7 drop?

table_1_down_ODC <- table_1_down  |> dplyr::filter(celltype == "ODC") #280
table_1_down_MG <- table_1_down  |> dplyr::filter(celltype == "MG") #222
table_1_down_OPC <- table_1_down  |> dplyr::filter(celltype == "OPC") #133
table_1_down_INH <- table_1_down  |> dplyr::filter(celltype == "INH") #181
table_1_down_EX <- table_1_down  |> dplyr::filter(celltype == "EX") #336

table_1_down_ASC <- table_1_down  |> dplyr::filter(celltype == "ASC") #437
table_1_down_PER.END <- table_1_down  |> dplyr::filter(celltype == "PER.END") #1 drop?


#### Create gene ensembl sets ####
df_list<- list(table_1_up_ODC, table_1_up_MG ,table_1_up_OPC , table_1_up_INH,table_1_up_EX,
               table_1_up_ASC, table_1_up_PER.END,
               table_1_down_ODC, table_1_down_MG ,table_1_down_OPC , table_1_down_INH,table_1_down_EX,
               table_1_down_ASC, table_1_down_PER.END )

res_1 <- purrr::map(df_list, get_ensembl, gene, "gene")


morabito_geneList <- list(
    table_1_up_ODC = res_1[[1]]$gene_ensembl_id,
    table_1_up_MG = res_1[[2]]$gene_ensembl_id,
    table_1_up_OPC = res_1[[3]]$gene_ensembl_id,
    table_1_up_INH = res_1[[4]]$gene_ensembl_id,
    table_1_up_EX = res_1[[5]]$gene_ensembl_id,
    table_1_up_ASC= res_1[[6]]$gene_ensembl_id,
    table_1_up_PER.END= res_1[[7]]$gene_ensembl_id,

    table_1_down_ODC = res_1[[8]]$gene_ensembl_id,
    table_1_down_MG = res_1[[9]]$gene_ensembl_id,
    table_1_down_OPC = res_1[[10]]$gene_ensembl_id,
    table_1_down_INH = res_1[[11]]$gene_ensembl_id,
    table_1_down_EX = res_1[[12]]$gene_ensembl_id,
    table_1_down_ASC= res_1[[13]]$gene_ensembl_id,
    table_1_down_PER.END= res_1[[14]]$gene_ensembl_id

)


#### Perform enrichment analysis ####
morabito_enrichment <- gene_set_enrichment(
    morabito_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment")

morabito_depleted <- gene_set_enrichment(
    morabito_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = TRUE)

# 15 0.000000 1.00000000       Ab+       table_1_up_ODC enrichment     0.1
# 16 0.000000 1.00000000       Ab+        table_1_up_MG enrichment     0.1
# 17 0.000000 1.00000000       Ab+       table_1_up_OPC enrichment     0.1
# 18 0.000000 1.00000000       Ab+       table_1_up_INH enrichment     0.1
# 19 0.000000 1.00000000       Ab+        table_1_up_EX enrichment     0.1
# 20 0.000000 1.00000000       Ab+       table_1_up_ASC enrichment     0.1
# 21 0.000000 1.00000000       Ab+   table_1_up_PER.END enrichment     0.1
# 22 0.000000 1.00000000       Ab+     table_1_down_ODC enrichment     0.1
# 23 0.000000 1.00000000       Ab+      table_1_down_MG enrichment     0.1
# 24 0.000000 1.00000000       Ab+     table_1_down_OPC enrichment     0.1
# 25 0.000000 1.00000000       Ab+     table_1_down_INH enrichment     0.1
# 26 0.000000 1.00000000       Ab+      table_1_down_EX enrichment     0.1
# 27 4.600446 0.21503694       Ab+     table_1_down_ASC enrichment     0.1
# 28 0.000000 1.00000000       Ab+ table_1_down_PER.END enrichment     0.1
# 29 2.399519 0.14136566  next_Ab+       table_1_up_ODC enrichment     0.1
# 30 2.441059 0.20703905  next_Ab+        table_1_up_MG enrichment     0.1
# 31 3.016744 0.15077258  next_Ab+       table_1_up_OPC enrichment     0.1
# 32 2.484643 0.33876672  next_Ab+       table_1_up_INH enrichment     0.1
# 33 0.000000 1.00000000  next_Ab+        table_1_up_EX enrichment     0.1
# 34 2.010391 0.20068352  next_Ab+       table_1_up_ASC enrichment     0.1
# 35 0.000000 1.00000000  next_Ab+   table_1_up_PER.END enrichment     0.1
# 36 2.073307 0.26071960  next_Ab+     table_1_down_ODC enrichment     0.1
# 37 1.609382 0.46953689  next_Ab+      table_1_down_MG enrichment     0.1
# 38 1.865692 0.42200935  next_Ab+     table_1_down_OPC enrichment     0.1
# 39 0.000000 1.00000000  next_Ab+     table_1_down_INH enrichment     0.1
# 40 1.297147 0.66934640  next_Ab+      table_1_down_EX enrichment     0.1
# 41 3.160512 0.04568396  next_Ab+     table_1_down_ASC enrichment     0.1
# 42 0.000000 1.00000000  next_Ab+ table_1_down_PER.END enrichment     0.1



table_5  <- read_excel("raw-data/GeneSets/3_snATAC-seq/Table S5_snATACseq.xlsx",
                       col_names = TRUE, skip = 2)

##### enrichment plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/06_morabito_enriched.pdf"), width = 12)
gene_set_enrichment_plot(
    morabito_enrichment ,
    xlabs = unique(morabito_enrichment$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(morabito_enrichment$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
                                                                             "YlOrRd")))(50)),
    cex = 1.2
)

dev.off()


pdf(paste0(output_dir, "/06_morabito_depleted.pdf"), width = 12)
gene_set_enrichment_plot(
    morabito_depleted ,
    xlabs = unique(morabito_depleted$ID),
    PThresh = 12,
    ORcut = 1.30103,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(morabito_depleted$test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
                                                                             "YlOrRd")))(50)),
    cex = 1.2
)

dev.off()
###snRNAseq

#Table S1_snRNA-seq -> on the "Supplementary Data 1e" sheet.
#ODC for oligodendrocytes; EX for excitatory neurons;
#MG for microglia; ASC for astrocytes; INH for inhibitory neurons;
#OPC for Oligodendrocyte precursor cells
