#### load relevant packages ####
library('readxl')
library('spatialLIBD')
library('dplyr')
library('sessioninfo')
library('here')
library('scran')
library('purrr')

### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensembl_function.R'))

### load modeling results
load(here('processed-data','11_grey_matter_only','wholegenome',
          'Visium_IF_AD_modeling_results.Rdata'))

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

input_dir <- "raw-data/GeneSets/2_snRNA-seq/1_Mathys et al_PFC/Mathys et al.xlsx"
mathys_ex <- read_excel(input_dir, sheet = 'Ex', skip = 1)
mathys_in <- read_excel(input_dir, sheet = 'In', skip = 1)
mathys_ast <- read_excel(input_dir, sheet = 'Ast', skip = 1)
mathys_oli <- read_excel(input_dir, sheet = 'Oli', skip = 1)
mathys_opc <- read_excel(input_dir, sheet = 'Opc', skip = 1)
mathys_mic <- read_excel(input_dir, sheet = 'Mic', skip = 1)


change_col_names <- function(df){
    colnames(df)[1] <- 'gene_set_np_v_p'
    colnames(df)[12] <- 'gene_set_np_v_ep'
    colnames(df)[23] <- 'gene_set_ep_v_lp'
    df

}


df.list <- list(mathys_ex, mathys_in, mathys_ast, mathys_oli, mathys_opc, mathys_mic)
df.list <- purrr::map(df.list, change_col_names)


# colnames(df.list[[1]])
# [1] "gene_set_np_v_p"           "IndModel.adj.pvals...2"
# [3] "no.pathology.mean...3"     "pathology.mean"
# [5] "IndModel.FC...5"           "MixedModel.z...6"
# [7] "MixedModel.p...7"          "DEGs.Ind.Model...8"
# [9] "DEGs.Ind.Mix.models...9"   "...10"
# [11] "...11"                     "gene_set_np_v_ep"
# [13] "IndModel.adj.pvals...13"   "no.pathology.mean...14"
# [15] "early.pathology.mean...15" "IndModel.FC...16"
# [17] "MixedModel.z...17"         "MixedModel.p...18"
# [19] "DEGs.Ind.Model...19"       "DEGs.Ind.Mix.models...20"
# [21] "...21"                     "...22"
# [23] "gene_set_ep_v_lp"          "IndModel.adj.pvals...24"
# [25] "late.pathology.mean"       "early.pathology.mean...26"
# [27] "IndModel.FC...27"          "MixedModel.z...28"
# [29] "MixedModel.p...29"         "DEGs.Ind.Model...30"
# [31] "DEGs.Ind.Mix.models...31"


res_1 <- purrr::map(df.list, get_ensembl, gene_set_np_v_p, 'gene_set_np_v_p')
res_2 <- purrr::map(df.list, get_ensembl, gene_set_np_v_ep, 'gene_set_np_v_ep')
res_3 <- purrr::map(df.list, get_ensembl, gene_set_ep_v_lp, 'gene_set_ep_v_lp')

#mathys_ex, mathys_in, mathys_ast, mathys_oli, mathys_opc, mathys_mic

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

##NAs in ensembl IDs?

mathys_enrichment <- gene_set_enrichment(
    mathys_geneList,
    fdr_cut = 0.1,
    modeling_results = modeling_results,
    model_type = "enrichment")


#### results ####
#            OR      Pval      test          ID model_type fdr_cut
# 1   0.0000000 1.0000000      none   np_v_p_ex enrichment     0.1
# 2   0.0000000 1.0000000      none   np_v_p_in enrichment     0.1
# 3   0.0000000 1.0000000      none  np_v_p_ast enrichment     0.1
# 4   0.0000000 1.0000000      none  np_v_p_oli enrichment     0.1
# 5   0.0000000 1.0000000      none  np_v_p_opc enrichment     0.1
# 6   0.0000000 1.0000000      none  np_v_p_mic enrichment     0.1
# 7   0.0000000 1.0000000      none  np_v_ep_ex enrichment     0.1
# 8   0.0000000 1.0000000      none  np_v_ep_in enrichment     0.1
# 9   0.0000000 1.0000000      none np_v_ep_ast enrichment     0.1
# 10  0.0000000 1.0000000      none np_v_ep_oli enrichment     0.1
# 11  0.0000000 1.0000000      none np_v_ep_opc enrichment     0.1
# 12  0.0000000 1.0000000      none np_v_ep_mic enrichment     0.1
# 13  0.0000000 1.0000000      none  ep_v_lp_ex enrichment     0.1
# 14  0.0000000 1.0000000      none  ep_v_lp_in enrichment     0.1
# 15  0.0000000 1.0000000      none ep_v_lp_ast enrichment     0.1
# 16  0.0000000 1.0000000      none ep_v_lp_oli enrichment     0.1
# 17  0.0000000 1.0000000      none ep_v_lp_opc enrichment     0.1
# 18  0.0000000 1.0000000      none ep_v_lp_mic enrichment     0.1
# 19  0.4788893 0.2967148       Ab+   np_v_p_ex enrichment     0.1
# 20  0.5183999 0.3272802       Ab+   np_v_p_in enrichment     0.1
# 21  1.3932396 1.0000000       Ab+  np_v_p_ast enrichment     0.1
# 22  0.5420900 0.3452591       Ab+  np_v_p_oli enrichment     0.1
# 23  1.4630534 1.0000000       Ab+  np_v_p_opc enrichment     0.1
# 24  1.2260288 1.0000000       Ab+  np_v_p_mic enrichment     0.1
# 25  0.4719318 0.2912762       Ab+  np_v_ep_ex enrichment     0.1
# 26  0.5294179 0.3356477       Ab+  np_v_ep_in enrichment     0.1
# 27  0.6245373 0.6335279       Ab+ np_v_ep_ast enrichment     0.1
# 28  0.5559702 0.3556873       Ab+ np_v_ep_oli enrichment     0.1
# 29  0.6629211 0.6421341       Ab+ np_v_ep_opc enrichment     0.1
# 30  1.4678441 1.0000000       Ab+ np_v_ep_mic enrichment     0.1
# 31  0.4917962 0.3067621       Ab+  ep_v_lp_ex enrichment     0.1
# 32  0.5526294 0.3531866       Ab+  ep_v_lp_in enrichment     0.1
# 33  0.7315388 0.6591440       Ab+ ep_v_lp_ast enrichment     0.1
# 34  1.4194631 1.0000000       Ab+ ep_v_lp_oli enrichment     0.1
# 35  0.7769018 0.6711408       Ab+ ep_v_lp_opc enrichment     0.1
# 36  2.0485276 0.4996621       Ab+ ep_v_lp_mic enrichment     0.1
# 37  1.2615313 0.8288620  next_Ab+   np_v_p_ex enrichment     0.1
# 38  1.1124907 1.0000000  next_Ab+   np_v_p_in enrichment     0.1
# 39  0.9354294 0.8431778  next_Ab+  np_v_p_ast enrichment     0.1
# 40  1.1634070 1.0000000  next_Ab+  np_v_p_oli enrichment     0.1
# 41  0.9825765 1.0000000  next_Ab+  np_v_p_opc enrichment     0.1
# 42  0.7648417 0.4223911  next_Ab+  np_v_p_mic enrichment     0.1
# 43  0.8474671 0.6631201  next_Ab+  np_v_ep_ex enrichment     0.1
# 44  0.9511741 0.8359926  next_Ab+  np_v_ep_in enrichment     0.1
# 45  1.1229084 1.0000000  next_Ab+ np_v_ep_ast enrichment     0.1
# 46  1.4652813 0.5412362  next_Ab+ np_v_ep_oli enrichment     0.1
# 47  0.8836409 0.7017850  next_Ab+ np_v_ep_opc enrichment     0.1
# 48  0.7045134 0.2217066  next_Ab+ np_v_ep_mic enrichment     0.1
# 49  1.0551167 1.0000000  next_Ab+  ep_v_lp_ex enrichment     0.1
# 50  0.9931209 1.0000000  next_Ab+  ep_v_lp_in enrichment     0.1
# 51  1.1245464 1.0000000  next_Ab+ ep_v_lp_ast enrichment     0.1
# 52  0.7257262 0.3302513  next_Ab+ ep_v_lp_oli enrichment     0.1
# 53  1.0363610 1.0000000  next_Ab+ ep_v_lp_opc enrichment     0.1
# 54  0.9061305 0.7715908  next_Ab+ ep_v_lp_mic enrichment     0.1
# 55  0.0000000 1.0000000       pT+   np_v_p_ex enrichment     0.1
# 56  0.0000000 1.0000000       pT+   np_v_p_in enrichment     0.1
# 57  0.0000000 1.0000000       pT+  np_v_p_ast enrichment     0.1
# 58  0.0000000 1.0000000       pT+  np_v_p_oli enrichment     0.1
# 59  0.0000000 1.0000000       pT+  np_v_p_opc enrichment     0.1
# 60  0.0000000 1.0000000       pT+  np_v_p_mic enrichment     0.1
# 61  0.0000000 1.0000000       pT+  np_v_ep_ex enrichment     0.1
# 62  0.0000000 1.0000000       pT+  np_v_ep_in enrichment     0.1
# 63  0.0000000 1.0000000       pT+ np_v_ep_ast enrichment     0.1
# 64  0.0000000 1.0000000       pT+ np_v_ep_oli enrichment     0.1
# 65  0.0000000 1.0000000       pT+ np_v_ep_opc enrichment     0.1
# 66  0.0000000 1.0000000       pT+ np_v_ep_mic enrichment     0.1
# 67  0.0000000 1.0000000       pT+  ep_v_lp_ex enrichment     0.1
# 68  0.0000000 1.0000000       pT+  ep_v_lp_in enrichment     0.1
# 69  0.0000000 1.0000000       pT+ ep_v_lp_ast enrichment     0.1
# 70  0.0000000 1.0000000       pT+ ep_v_lp_oli enrichment     0.1
# 71  0.0000000 1.0000000       pT+ ep_v_lp_opc enrichment     0.1
# 72  0.0000000 1.0000000       pT+ ep_v_lp_mic enrichment     0.1
# 73  0.0000000 1.0000000  next_pT+   np_v_p_ex enrichment     0.1
# 74  0.0000000 1.0000000  next_pT+   np_v_p_in enrichment     0.1
# 75  0.0000000 1.0000000  next_pT+  np_v_p_ast enrichment     0.1
# 76  0.0000000 1.0000000  next_pT+  np_v_p_oli enrichment     0.1
# 77  0.0000000 1.0000000  next_pT+  np_v_p_opc enrichment     0.1
# 78  0.0000000 1.0000000  next_pT+  np_v_p_mic enrichment     0.1
# 79  0.0000000 1.0000000  next_pT+  np_v_ep_ex enrichment     0.1
# 80  0.0000000 1.0000000  next_pT+  np_v_ep_in enrichment     0.1
# 81  0.0000000 1.0000000  next_pT+ np_v_ep_ast enrichment     0.1
# 82  0.0000000 1.0000000  next_pT+ np_v_ep_oli enrichment     0.1
# 83  0.0000000 1.0000000  next_pT+ np_v_ep_opc enrichment     0.1
# 84  0.0000000 1.0000000  next_pT+ np_v_ep_mic enrichment     0.1
# 85  0.0000000 1.0000000  next_pT+  ep_v_lp_ex enrichment     0.1
# 86  0.0000000 1.0000000  next_pT+  ep_v_lp_in enrichment     0.1
# 87  0.0000000 1.0000000  next_pT+ ep_v_lp_ast enrichment     0.1
# 88  0.0000000 1.0000000  next_pT+ ep_v_lp_oli enrichment     0.1
# 89  0.0000000 1.0000000  next_pT+ ep_v_lp_opc enrichment     0.1
# 90  0.0000000 1.0000000  next_pT+ ep_v_lp_mic enrichment     0.1
# 91  0.0000000 1.0000000      both   np_v_p_ex enrichment     0.1
# 92  0.0000000 1.0000000      both   np_v_p_in enrichment     0.1
# 93  0.0000000 1.0000000      both  np_v_p_ast enrichment     0.1
# 94  0.0000000 1.0000000      both  np_v_p_oli enrichment     0.1
# 95  0.0000000 1.0000000      both  np_v_p_opc enrichment     0.1
# 96  0.0000000 1.0000000      both  np_v_p_mic enrichment     0.1
# 97  0.0000000 1.0000000      both  np_v_ep_ex enrichment     0.1
# 98  0.0000000 1.0000000      both  np_v_ep_in enrichment     0.1
# 99  0.0000000 1.0000000      both np_v_ep_ast enrichment     0.1
# 100 0.0000000 1.0000000      both np_v_ep_oli enrichment     0.1
# 101 0.0000000 1.0000000      both np_v_ep_opc enrichment     0.1
# 102 0.0000000 1.0000000      both np_v_ep_mic enrichment     0.1
# 103 0.0000000 1.0000000      both  ep_v_lp_ex enrichment     0.1
# 104 0.0000000 1.0000000      both  ep_v_lp_in enrichment     0.1
# 105 0.0000000 1.0000000      both ep_v_lp_ast enrichment     0.1
# 106 0.0000000 1.0000000      both ep_v_lp_oli enrichment     0.1
# 107 0.0000000 1.0000000      both ep_v_lp_opc enrichment     0.1
# 108 0.0000000 1.0000000      both ep_v_lp_mic enrichment     0.1
# 109 0.0000000 1.0000000 next_both   np_v_p_ex enrichment     0.1
# 110 0.0000000 1.0000000 next_both   np_v_p_in enrichment     0.1
# 111 0.0000000 1.0000000 next_both  np_v_p_ast enrichment     0.1
# 112 0.0000000 1.0000000 next_both  np_v_p_oli enrichment     0.1
# 113 0.0000000 1.0000000 next_both  np_v_p_opc enrichment     0.1
# 114 0.0000000 1.0000000 next_both  np_v_p_mic enrichment     0.1
# 115 0.0000000 1.0000000 next_both  np_v_ep_ex enrichment     0.1
# 116 0.0000000 1.0000000 next_both  np_v_ep_in enrichment     0.1
# 117 0.0000000 1.0000000 next_both np_v_ep_ast enrichment     0.1
# 118 0.0000000 1.0000000 next_both np_v_ep_oli enrichment     0.1
# 119 0.0000000 1.0000000 next_both np_v_ep_opc enrichment     0.1
# 120 0.0000000 1.0000000 next_both np_v_ep_mic enrichment     0.1
# 121 0.0000000 1.0000000 next_both  ep_v_lp_ex enrichment     0.1
# 122 0.0000000 1.0000000 next_both  ep_v_lp_in enrichment     0.1
# 123 0.0000000 1.0000000 next_both ep_v_lp_ast enrichment     0.1
# 124 0.0000000 1.0000000 next_both ep_v_lp_oli enrichment     0.1
# 125 0.0000000 1.0000000 next_both ep_v_lp_opc enrichment     0.1
# 126 0.0000000 1.0000000 next_both ep_v_lp_mic enrichment     0.1

##### enrichment plotting #####
output_dir <- here("plots", "14_external_gene_sets")
pdf(paste0(output_dir, "/03_mathys.pdf"), width = 15)
gene_set_enrichment_plot(
    mathys_enrichment,
    xlabs = unique(mathys_enrichment$ID),
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    layerHeights = c(0, seq_len(length(unique(mathys_enrichment $test)))) * 15,
    mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
                                                                             "YlOrRd")))(50)),
    cex = 1.2
)

dev.off()



