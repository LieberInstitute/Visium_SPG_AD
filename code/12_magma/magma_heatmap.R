# library(sgejobs)

# sgejobs::job_single(
#     "magma_heatmap",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript magma_heatmap.R",
#     create_logdir = TRUE,
#     task_num = 3
# )

# library(rtracklayer)
# library(GenomicRanges)
# library(jaffelab)
# library(SingleCellExperiment)
# library(fields)
library(reshape2)
library(RColorBrewer)
library(grid)
library(sessioninfo)
library(here)
library(dplyr)
here()

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")

magmaStats <- list()

gene_num = c("50",  "100", "200")


magmaStats[["ITC"]][["AD.Jansen.2019"]] <- read.table(here("code", "12_magma", "01_Jansen_2019", paste0("ad_gwas_", gene_num[SGE_TASK_ID], ".gsa.out")), header = T)
magmaStats[["ITC"]][["FTD.Ferrari.2014"]] <- read.table(here("code", "12_magma", "02_Lancet_2014", paste0("ftd_gwas_", gene_num[SGE_TASK_ID], ".gsa.out")), header = T)
magmaStats[["ITC"]][["PD.Nalls.2019"]] <- read.table(here("code", "12_magma", "03_Nalls_2019", paste0("pd_gwas_", gene_num[SGE_TASK_ID], ".gsa.out")), header = T)


## Merge to assess significance thresholds ===
magmaStats_list <- lapply(magmaStats, function(m) {
    z <- sapply(m, "[[", "P")
    rownames(z) <- m[[1]]$VARIABLE
    z
})

magmaStats_wide <- as.data.frame(do.call("rbind", magmaStats_list))

magmaStats_wide$Region <- rep("ITC",
    times = sapply(magmaStats_list, nrow)
)
magmaStats_wide$Pathology_Type <- rownames(magmaStats_wide)


##Compute disorder-wise FDRs and Bonferronis
magmaStats_wide$FDR_AD <- p.adjust(magmaStats_wide$AD.Jansen.2019, "fdr")
magmaStats_wide$FDR_FTD <- p.adjust(magmaStats_wide$FTD.Ferrari.2014, "fdr")
magmaStats_wide$FDR_PD <- p.adjust(magmaStats_wide$PD.Nalls.2019, "fdr")

magmaStats_wide$Bonf_AD <- p.adjust(magmaStats_wide$AD.Jansen.2019, "bonf")
magmaStats_wide$Bonf_FTD <- p.adjust(magmaStats_wide$FTD.Ferrari.2014, "bonf")
magmaStats_wide$Bonf_PD <- p.adjust(magmaStats_wide$PD.Nalls.2019, "bonf")

#### TOP 50 GENES ####

#           AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019 Region Pathology_Type
# Ab               0.48500         0.793820       0.41681    ITC             Ab
# both             0.73405         0.040561       0.17218    ITC           both
# next_Ab          0.98361         0.620600       0.48718    ITC        next_Ab
# next_both        0.21078         0.610990       0.26610    ITC      next_both
# next_pT          0.58011         0.054698       0.88031    ITC        next_pT
# none             0.11293         0.254650       0.26922    ITC           none
# pT               0.64360         0.352770       0.40913    ITC             pT
#               FDR_AD   FDR_FTD    FDR_PD Bonf_AD Bonf_FTD Bonf_PD
# Ab        0.8563917 0.7938200 0.5683767 1.00000 1.000000       1
# both      0.8563917 0.1914430 0.5683767 1.00000 0.283927       1
# next_Ab   0.9836100 0.7240333 0.5683767 1.00000 1.000000       1
# next_both 0.7377300 0.7240333 0.5683767 1.00000 1.000000       1
# next_pT   0.8563917 0.1914430 0.8803100 1.00000 0.382886       1
# none      0.7377300 0.5941833 0.5683767 0.79051 1.000000       1
# pT        0.8563917 0.6173475 0.5683767 1.00000 1.000000       1


#### TOP 100 GENES ####
#               AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019 Region Pathology_Type
# Ab              0.289260         0.765570      0.272020    ITC             Ab
# both            0.696740         0.010701      0.067224    ITC           both
# next_Ab         0.950070         0.677650      0.245110    ITC        next_Ab
# next_both       0.288460         0.791890      0.278450    ITC      next_both
# next_pT         0.691490         0.139700      0.645650    ITC        next_pT
# none            0.014898         0.289230      0.031440    ITC           none
# pT              0.445470         0.746460      0.555750    ITC             pT
#           FDR_AD  FDR_FTD   FDR_PD  Bonf_AD Bonf_FTD  Bonf_PD
# Ab        0.6749400 0.791890 0.389830 1.000000 1.000000 1.000000
# both      0.8128633 0.074907 0.235284 1.000000 0.074907 0.470568
# next_Ab   0.9500700 0.791890 0.389830 1.000000 1.000000 1.000000
# next_both 0.6749400 0.791890 0.389830 1.000000 1.000000 1.000000
# next_pT   0.8128633 0.488950 0.645650 1.000000 0.977900 1.000000
# none      0.1042860 0.674870 0.220080 0.104286 1.000000 0.220080
# pT        0.7795725 0.791890 0.645650 1.000000 1.000000 1.000000

#### TOP 200 GENES ####
#               AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019 Region Pathology_Type
# Ab               0.16333         0.692780      0.182410    ITC             Ab
# both             0.87265         0.091237      0.052265    ITC           both
# next_Ab          0.87754         0.113850      0.297910    ITC        next_Ab
# next_both        0.18025         0.477900      0.188100    ITC      next_both
# next_pT          0.72660         0.219960      0.331360    ITC        next_pT
# none             0.21561         0.048492      0.080990    ITC           none
# pT               0.53431         0.521780      0.638230    ITC             pT
#           FDR_AD   FDR_FTD    FDR_PD Bonf_AD Bonf_FTD  Bonf_PD
# Ab        0.50309 0.6927800 0.3291750       1 1.000000 1.000000
# both      0.87754 0.2656500 0.2834650       1 0.638659 0.365855
# next_Ab   0.87754 0.2656500 0.3865867       1 0.796950 1.000000
# next_both 0.50309 0.6087433 0.3291750       1 1.000000 1.000000
# next_pT   0.87754 0.3849300 0.3865867       1 1.000000 1.000000
# none      0.50309 0.2656500 0.2834650       1 0.339444 0.566930
# pT        0.87754 0.6087433 0.6382300       1 1.000000 1.000000



#### reshape to long ####
magmaStats_long <- reshape2::melt(magmaStats_wide[, 1:5])
colnames(magmaStats_long)[3:4] <- c("GWAS", "P")
dim(magmaStats_long)

## Calculate Global FDR and Bonf, i.e. across all Pathologies X Disorders ##
magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
magmaStats_long$P.adj.bonf <- p.adjust(magmaStats_long$P, "bonf")


#### TOP 200 GENES ####
# Region Pathology_Type             GWAS        P P.adj.fdr P.adj.bonf
# 1     ITC             Ab   AD.Jansen.2019 0.163330 0.4199236          1
# 2     ITC           both   AD.Jansen.2019 0.872650 0.8775400          1
# 3     ITC        next_Ab   AD.Jansen.2019 0.877540 0.8775400          1
# 4     ITC      next_both   AD.Jansen.2019 0.180250 0.4199236          1
# 5     ITC        next_pT   AD.Jansen.2019 0.726600 0.8030842          1
# 6     ITC           none   AD.Jansen.2019 0.215610 0.4199236          1
# 7     ITC             pT   AD.Jansen.2019 0.534310 0.7012819          1
# 8     ITC             Ab FTD.Ferrari.2014 0.692780 0.8030842          1
# 9     ITC           both FTD.Ferrari.2014 0.091237 0.4199236          1
# 10    ITC        next_Ab FTD.Ferrari.2014 0.113850 0.4199236          1
# 11    ITC      next_both FTD.Ferrari.2014 0.477900 0.7012819          1
# 12    ITC        next_pT FTD.Ferrari.2014 0.219960 0.4199236          1
# 13    ITC           none FTD.Ferrari.2014 0.048492 0.4199236          1
# 14    ITC             pT FTD.Ferrari.2014 0.521780 0.7012819          1
# 15    ITC             Ab    PD.Nalls.2019 0.182410 0.4199236          1
# 16    ITC           both    PD.Nalls.2019 0.052265 0.4199236          1
# 17    ITC        next_Ab    PD.Nalls.2019 0.297910 0.5213425          1
# 18    ITC      next_both    PD.Nalls.2019 0.188100 0.4199236          1
# 19    ITC        next_pT    PD.Nalls.2019 0.331360 0.5352738          1
# 20    ITC           none    PD.Nalls.2019 0.080990 0.4199236          1
# 21    ITC             pT    PD.Nalls.2019 0.638230 0.7884018          1

#### Extract and shape all the Beta values from the magma output ####
magmaStats_list.beta <- lapply(magmaStats, function(m) {
    z <- sapply(m, "[[", "BETA")
    rownames(z) <- m[[1]]$VARIABLE
    z
})
magmaStats_wide.beta <- as.data.frame(do.call("rbind", magmaStats_list.beta))
magmaStats_wide.beta$Region <- rep("ITC",
    times = sapply(magmaStats_list.beta, nrow)
)
magmaStats_wide.beta$Pathology_Type <- rownames(magmaStats_wide.beta)

####betas for top 200 ####
# AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019 Region Pathology_Type
# Ab             0.0542280       -0.0384820      0.055908    ITC             Ab
# both          -0.0661240        0.1026800      0.106150    ITC           both
# next_Ab       -0.0672980        0.0895870      0.032954    ITC        next_Ab
# next_both      0.0502020        0.0038814      0.054061    ITC      next_both
# next_pT       -0.0345860        0.0554810      0.028454    ITC        next_pT
# none           0.0432630        0.1114300      0.086360    ITC           none
# pT            -0.0049454       -0.0037398     -0.022830    ITC             pT


magmaStats_long.beta <- reshape2::melt(magmaStats_wide.beta)
colnames(magmaStats_long.beta)[3:4] <- c("GWAS", "Beta")

# Check before appending
table(paste0(magmaStats_long$Pathology_Type, ":", magmaStats_long$GWAS) ==
    paste0(magmaStats_long.beta$Pathology_Type, ":", magmaStats_long.beta$GWAS))

# cbind Beta
magmaStats_long$Beta <- magmaStats_long.beta$Beta
# Reorder
magmaStats_long <- magmaStats_long[, c("Region", "Pathology_Type", "GWAS", "Beta", "P", "P.adj.fdr")]

magmaStats_wide$AD_Beta <- magmaStats_wide.beta$AD.Jansen.2019
magmaStats_wide$FTD_Beta <- magmaStats_wide.beta$FTD.Ferrari.2014
magmaStats_wide$PD_Beta <- magmaStats_wide.beta$PD.Nalls.2019


## Print to CSV ===
write.csv(magmaStats_long,
    file = here(
        "code", "12_magma",
        "Global_FDRs_Bonf_3GWAS_x_7_Pathologies.csv"
    ),
    row.names = F, quote = F
)

write.csv(magmaStats_wide,
          file = here(
              "code", "12_magma",
              "GWAS_wise_FDRs_Bonf_3GWAS_x_7_Pathologies.csv"
          ),
          row.names = F, quote = F
)

head(magmaStats_long)

#### Make heatmap ####
midpoint <- function(x) x[-length(x)] + diff(x) / 2

MAGMAplot <- function(region, Pthresh, fdrThresh, ...) {
    ####what does magmaStats[[region]] look like? ####
    # $AD.Jansen.2019
    # VARIABLE TYPE NGENES       BETA    BETA_STD       SE       P
    # 1        Ab  SET    194  0.0542280  0.00428220 0.055285 0.16333
    # 2      both  SET    190 -0.0661240 -0.00516780 0.058053 0.87265
    # 3   next_Ab  SET    188 -0.0672980 -0.00523200 0.057875 0.87754
    # 4 next_both  SET    196  0.0502020  0.00398460 0.054901 0.18025
    # 5   next_pT  SET    193 -0.0345860 -0.00272410 0.057396 0.72660
    # 6      none  SET    190  0.0432630  0.00338110 0.054963 0.21561
    # 7        pT  SET    193 -0.0049454 -0.00038952 0.057434 0.53431
    #
    # $FTD.Ferrari.2014
    # VARIABLE TYPE NGENES       BETA    BETA_STD       SE        P
    # 1        Ab  SET    192 -0.0384820 -0.00306540 0.076391 0.692780
    # 2      both  SET    187  0.1026800  0.00807250 0.077015 0.091237
    # 3   next_Ab  SET    187  0.0895870  0.00704340 0.074263 0.113850
    # 4 next_both  SET    195  0.0038814  0.00031157 0.070018 0.477900
    # 5   next_pT  SET    190  0.0554810  0.00439660 0.071834 0.219960
    # 6      none  SET    187  0.1114300  0.00876030 0.067135 0.048492
    # 7        pT  SET    192 -0.0037398 -0.00029790 0.068475 0.521780
    #
    # $PD.Nalls.2019
    # VARIABLE TYPE NGENES      BETA   BETA_STD       SE        P
    # 1        Ab  SET    194  0.055908  0.0044136 0.061692 0.182410
    # 2      both  SET    190  0.106150  0.0082937 0.065391 0.052265
    # 3   next_Ab  SET    188  0.032954  0.0025612 0.062126 0.297910
    # 4 next_both  SET    196  0.054061  0.0042896 0.061091 0.188100
    # 5   next_pT  SET    192  0.028454  0.0022347 0.065237 0.331360
    # 6      none  SET    190  0.086360  0.0067473 0.061752 0.080990
    # 7        pT  SET    191 -0.022830 -0.0017884 0.064542 0.638230

    #### Set up -log10(p_value) ####
    wide_p <- sapply(magmaStats[[region]], function(x) {
        cbind(-log10(x$P))
    })
    rownames(wide_p) <- magmaStats[[region]][[1]]$VARIABLE
    wide_p[wide_p > Pthresh] <- Pthresh
    wide_p <- round(wide_p[rev(sort(rownames(wide_p))), ], 3)

    #### wide_p for top 200 genes####
    # AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019
    # pT                 0.272            0.283         0.195
    # none               0.666            1.314         1.092
    # next_pT            0.139            0.658         0.480
    # next_both          0.744            0.321         0.726
    # next_Ab            0.057            0.944         0.526
    # both               0.059            1.040         1.282
    # Ab                 0.787            0.159         0.739
    #### Set up betas ####
    wide_beta <- sapply(magmaStats[[region]], function(x) {
        cbind(x$BETA)
    })
    rownames(wide_beta) <- magmaStats[[region]][[1]]$VARIABLE
    wide_beta <- round(wide_beta[rev(sort(rownames(wide_beta))), ], 2)
    # beta = 0.11 for both
    ####what does wide_beta look like for top 200? ####
    # AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019
    # pT                  0.00             0.00         -0.02
    # none                0.04             0.11          0.09
    # next_pT            -0.03             0.06          0.03
    # next_both           0.05             0.00          0.05
    # next_Ab            -0.07             0.09          0.03
    # both               -0.07             0.10          0.11
    # Ab                  0.05            -0.04          0.06

    #### Use empirical cutoff (independent p=0.05) for printing betas ####
    wide_beta[wide_p < -log10(0.05)] <- ""
    # anything with independent p-value lower than 0.05 will not be printed in the
    # heatmap.

    # wide_beta
    # AD.Jansen.2019 FTD.Ferrari.2014 PD.Nalls.2019
    # pT        ""             ""               ""
    # none      ""             "0.11"           ""
    # next_pT   ""             ""               ""
    # next_both ""             ""               ""
    # next_Ab   ""             ""               ""
    # both      ""             ""               ""
    # Ab        ""             ""               ""

    # and Bonf cutoff for bolding
    customFont <- ifelse(wide_p < -log10(fdrThresh), 1, 2)
    customCex <- ifelse(wide_p < -log10(fdrThresh), 0.9, 1.0)


    ## Plot
    clusterHeights <- seq(0, 160, length.out = nrow(wide_p) + 1)
    mypal <- c("white", colorRampPalette(brewer.pal(9, "YlOrRd"))(60))[1:30]
    xlabs <- colnames(wide_p)

    # Heatmap of p's
    image.plot(
        x = seq(0, ncol(wide_p), by = 1), y = clusterHeights, z = as.matrix(t(wide_p)),
        col = mypal, xaxt = "n", yaxt = "n", xlab = "", ylab = ""
    )
    axis(2, rownames(wide_p), at = midpoint(clusterHeights), las = 1)
    axis(1, rep("", ncol(wide_p)), at = seq(0.5, ncol(wide_p) - 0.5))
    text(
        x = seq(0.5, ncol(wide_p) - 0.5), y = -1 * max(nchar(xlabs)) / 2, xlabs,
        xpd = TRUE, srt = 45, cex = 1.2, adj = 1
    )
    abline(h = clusterHeights, v = 0:ncol(wide_p))

    # Print top decile of betas
    text(
        x = rep(seq(0.5, ncol(wide_p) - 0.5), each = nrow(wide_p)),
        y = rep(midpoint(clusterHeights), ncol(wide_p)),
        as.character(wide_beta),
        # If Bonf, a little bigger
        cex = customCex,
        # If Bonf, 2 (bold)
        font = customFont
    )
}

# Plot
dir_output <- dir.create(here("plots", "magma"))
pdf(here("plots", "magma", paste0("heatmap_", "top_", gene_num[SGE_TASK_ID], "_genes", ".pdf"), w = 6))
par(mar = c(8.5, 7.5, 6, 1), cex.axis = 1.0, cex.lab = 0.5)
MAGMAplot(region = "ITC", Pthresh = 12, fdrThresh = 1) ## set to 1 because didn't want to make a cutoff based on FDR/Bonf
abline(v = 7, lwd = 3)
text(x = 1.5, y = 185, "MAGMA gene set analyses: ITC pathology types", xpd = TRUE, cex = 1.5, font = 2)
text(x = 1.5, y = 175, "(Betas for empirically significant associations, p < 0.05)", xpd = TRUE, cex = 1, font = 1)
grid::grid.text(label = "-log10(p-value)", x = 0.92, y = 0.825, gp = gpar(fontsize = 9))
dev.off()


#
# ##### Explore gene-level MAGMA stats =========
#
# spe_wholegenome <- readRDS("/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/08_harmony_BayesSpace/wholegenome/spe_harmony_wholegenome.rds")
# pathology_types <- read.csv(here(
#     "processed-data", "09_pathology_vs_BayesSpace",
#     "pathology_levels", "path_groups", "clusters.csv"
# ))
# colData(spe_wholegenome)$path_type <- pathology_types$cluster
#
# markers <- as_tibble(read.table(here("code", "magma", "pvalues_top_200.txt")))
#
# markers$gene_name <- rowData(spe_wholegenome)$gene_name[match(markers$V2, rowData(spe_wholegenome)$gene_id)]
#
# #### AD Genes
#
# ad.genes.out <- read.table("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/SNP_Data/AD_Jansen2019_LC_snp-wise.genes.out",
#     sep = "", header = T
# )
#
# head(ad.genes.out)
#
# round(quantile(ad.genes.out$ZSTAT, probs = seq(0.9, 1, by = 0.01)), 3)
# # 90%    91%    92%    93%    94%    95%    96%    97%    98%    99%   100%
# # 1.535  1.603  1.679  1.760  1.860  1.990  2.131  2.324  2.638  3.282 10.183
#
#
# round(quantile(ad.genes.out$ZSTAT, probs = seq(0.95, 1, by = 0.0025)), 3)
# # 95% 95.25%  95.5% 95.75%    96% 96.25%  96.5% 96.75%    97% 97.25%  97.5%
# # 1.990  2.020  2.054  2.091  2.131  2.166  2.219  2.267  2.324  2.384  2.444
# # 97.75%    98% 98.25%  98.5% 98.75%    99% 99.25%  99.5% 99.75%   100%
# #     2.542  2.638  2.734  2.850  3.038  3.282  3.606  4.049  4.908 10.183
#
# jansen_sig_loci <- c(
#     "ADAMTS4", "CR1", "BIN1", "INPPD5", "HESX1", "CLNK", "HLA-DRB1", "TREM2", "CD2AP",
#     "ZCWPW1", "EPHA1", "CNTNAP2", "CLU/PTK2B", "ECHDC3", "MS4A6A", "PICALM", "SORL1", "SLC24A4",
#     "ADAM10", "APH1B", "KAT8", "SCIMP", "ABI3", "ALPK2", "ABCA7", "APOE",
#     "AC074212.3", "CD33", "CASS4"
# )
# table(jansen_sig_loci %in% rowData(spe_wholegenome)$gene_name)
# # FALSE  TRUE
# # 3    26
# #######
#
# setdiff(jansen_sig_loci, rowData(spe_wholegenome)$gene_name)
# # "INPPD5"     "CLU/PTK2B"  "AC074212.3"
#
# genes_for_plotting <- intersect(jansen_sig_loci, rowData(spe_wholegenome)$gene_name)
#
#
#
# ## Colorblind colors
# colors_pathology <- setNames(
#     c(
#         "grey90",
#         paletteer::paletteer_d("dichromat::SteppedSequential_5")[rep(c(6, 18), each = 2) + c(0, 3)],
#         paletteer::paletteer_d("beyonce::X7")[4:5]
#     )[c(1:3, 6:7, 4:5)],
#     c("none", "Ab+", "next_Ab+", "pT+", "next_pT+", "both", "next_both")
# )
# Footer
#
# colors2print <- colors_pathology
#
#
# dir_output <- dir.create(here("plots", "violin_plots"))
# pdf(here("plots", "violin_plots", "MAGMA_reported_ad_loci.pdf"), width = 9, height = 8)
# plotExpressionCustom(
#     sce = spe_wholegenome,
#     exprs_values = "logcounts",
#     features = genes_for_plotting,
#     features_name = "",
#     anno_name = "path_type",
#     ncol = 4, point_alpha = 0.4,
#     scales = "free_y", swap_rownames = "gene_name"
# ) +
#     scale_color_manual(values = colors2print) +
#     ggtitle(label = paste0("AD risk loci expression across pathology domains in ITC")) +
#     theme(plot.title = element_text(size = 12))
# dev.off()
#
#
#
# round(quantile(pd.gene.out$ZSTAT, probs = seq(0.9, 1, by = 0.01)), 3)
# # 90%   91%   92%   93%   94%   95%   96%   97%   98%   99%  100%
# # 1.368 1.448 1.532 1.629 1.747 1.868 2.016 2.219 2.473 3.048 7.916
#
# round(quantile(pd.gene.out$ZSTAT, probs = seq(0.95, 1, by = 0.0025)), 3)
# # 95% 95.25%  95.5% 95.75%    96% 96.25%  96.5% 96.75%    97% 97.25%  97.5%
# # 1.868  1.898  1.934  1.974  2.016  2.063  2.121  2.169  2.219  2.274  2.337
# # 97.75%    98% 98.25%  98.5% 98.75%    99% 99.25%  99.5% 99.75%   100%
# #     2.399  2.473  2.566  2.684  2.867  3.048  3.342  3.734  5.088  7.916
#
# table(pd.gene.out$ZSTAT >= 5.0) # 83
# # 83 (90 reported index SNPs)
# topGenes.pd <- pd.gene.out$GENE[pd.gene.out$ZSTAT >= 3.7]
# table(topGenes.pd %in% rowData(spe_wholegenome)$gene_id) # all there
# # FALSE  TRUE
# # 16   143
# topGenes.pd <- rowData(spe_wholegenome)[match(topGenes.pd, rowData(spe_wholegenome)$gene_id)]
# topGenes.pd
# # [1] "ADAM15"     "EFNA4"      "DCST1-AS1"  "EFNA3"      "EFNA1"
# # [6] "SLC50A1"    "DPM3"       "TRIM46"     "ASH1L-IT1"  "VAMP4"
# # [11] "SLC45A3"    "NUCKS1"     NA           "RAB29"      "SLC41A1"
# # [16] "IFT172"     "GCKR"       "MAP4K4"     "TMEM163"    "CCNT2-AS1"
# # [21] NA           "CCNT2"      "STK39"      "DNASE1L3"   "ABHD6"
# # [26] "MED12L"     "P2RY14"     "P2RY12"     "AC092953.2" "DCUN1D1"
# # [31] NA           "MCCC1"      "CPLX1"      "AC139887.3" NA
# # [36] "GAK"        "TMEM175"    "DGKQ"       "IDUA"       "SLC26A1"
# # [41] "FGFRL1"     NA           "FBXL5"      "FAM200B"    "BST1"
# # [46] "SCARB2"     "FAM47E"     "AC034139.1" "STBD1"      "CCDC158"
# # [51] "AC093866.1" NA           "AC097478.2" "SNCA"       "SNCA-AS1"
# # [56] "MMRN1"      "ELOVL7"     "AC104113.1" "ERCC8"      NA
# # [61] "NDUFAF2"    "SMIM15-AS1" "SMIM15"     "ZNF391"     "AL021918.3"
# # [66] "KLHL7"      "KLHL7-DT"   "NUP42"      "AC005082.1" "AC005082.2"
# # [71] "GPNMB"      "SH3GL2"     NA           "IGSF9B"     "SLC2A13"
# # [76] "LINC02471"  "LRRK2"      "AC079630.1" "MUC19"      NA
# # [81] "KNTC1"      "AC026333.4" NA           "AL031600.2" "TELO2"
# # [86] "LINC02124"  "NDUFB10"    "MSRB1"      "PRR14"      "FBRS"
# # [91] "SRCAP"      "AC093249.6" "TMEM265"    "ZNF629"     "AC106886.3"
# # [96] "BCL7C"      "MIR762HG"   "CTF1"       "AC135048.2" "FBXL19"
# # [101] "FBXL19-AS1" "ORAI3"      "SETD1A"     "AC135048.3" "HSD3B7"
# # [106] "STX1B"      "STX4"       "AC135050.3" "AC135050.1" NA
# # [111] "ZNF646"     "ZNF668"     "BCKDK"      "PRSS53"     "KAT8"
# # [116] "VKORC1"     "AC135050.5" "AC135050.6" "PRSS8"      "TOX3"
# # [121] "CASC16"     "AC026462.3" "OR3A3"      "TTC19"      "ZSWIM7"
# # [126] "AC002553.1" "AC002553.2" "NCOR1"      "PIGL"       "AC003070.2"
# # [131] NA           "ARHGAP27"   NA           "PLEKHM1"    NA
# # [136] NA           "CRHR1"      "MAPT-AS1"   "SPPL2C"     "MAPT"
# # [141] "MAPT-IT1"   "CR936218.2" "STH"        "KANSL1"     "KANSL1-AS1"
# # [146] "LRRC37A"    "ARL17B"     "NSF"        "WNT3"       "RIT2"
# # [151] "AC100779.1" "FASTKD5"    "LZTS3"      "ITPA"       "DDRGK1"
# # [156] NA           "AP001437.2" "DYRK1A"     "AP001437.1"
#
#
# load("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/SCE/sce_updated_LC.rda", verbose = T)
