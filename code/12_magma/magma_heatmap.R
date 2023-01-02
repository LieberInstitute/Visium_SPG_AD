# library(sgejobs)

# sgejobs::job_single(
#     "magma_heatmap",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript magma_heatmap.R",
#     create_logdir = TRUE,
#     task_num = 4
# )

# library(rtracklayer)
# library(GenomicRanges)
# library(jaffelab)
# library(SingleCellExperiment)
# library(fields)
library(reshape2)
library(RColorBrewer)
library(rtracklayer)
library(GenomicRanges)
library(jaffelab)
library(SingleCellExperiment)
library(fields)
library(RColorBrewer)
library(dplyr)
library(grid)
library(sessioninfo)
library(here)

here()

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")

magmaStats <- list()

k <-  as.numeric(Sys.getenv("SGE_TASK_ID"))
print(k)
gene_set <- c("50", "100", "200", "fdr")
print(gene_set[k])


magmaStats[["ITC"]][["AD.Jansen.2019"]] <- read.table(here("code", "12_magma", "01_Jansen_2019", paste0("ad_gwas_", gene_set[k], ".gsa.out")), header = T)
magmaStats[["ITC"]][["FTD.Ferrari.2014"]] <- read.table(here("code", "12_magma", "02_Lancet_2014", paste0("ftd_gwas_", gene_set[k], ".gsa.out")), header = T)
magmaStats[["ITC"]][["PD.Nalls.2019"]] <- read.table(here("code", "12_magma", "03_Nalls_2019", paste0("pd_gwas_", gene_set[k], ".gsa.out")), header = T)


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


## Compute disorder-wise FDRs and Bonferronis
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

#### betas for top 200 ####
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
        "code", "12_magma", "heatmaps",
        paste0("Global_FDRs_Bonf_3GWAS_x_7_Pathologies_", gene_set[k],".csv")),
    row.names = F, quote = F
)

write.csv(magmaStats_wide,
    file = here(
        "code", "12_magma","heatmaps",
        paste0("GWAS_wise_FDRs_Bonf_3GWAS_x_7_Pathologies_", gene_set[k],".csv")),
    row.names = F, quote = F
)

head(magmaStats_long)

#### Make heatmap ####
midpoint <- function(x) x[-length(x)] + diff(x) / 2

MAGMAplot <- function(region, Pthresh, fdrThresh, ...) {
    #### what does magmaStats[[region]] look like? ####
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
    #### what does wide_beta look like for top 200? ####
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
    fields::image.plot(
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
pdf(here("plots", "magma", paste0("heatmap_", "top_", gene_set[k], "_genes", ".pdf")), w = 6)
par(mar = c(8.5, 7.5, 6, 1), cex.axis = 1.0, cex.lab = 0.5)
MAGMAplot(region = "ITC", Pthresh = 12, fdrThresh = 1) ## set to 1 because didn't want to make a cutoff based on FDR/Bonf
abline(v = 7, lwd = 3)
text(x = 1.5, y = 185, "MAGMA gene set analyses: ITC pathology types", xpd = TRUE, cex = 1.5, font = 2)
text(x = 1.5, y = 175, "(Betas for empirically significant associations, p < 0.05)", xpd = TRUE, cex = 1, font = 1)
grid::grid.text(label = "-log10(p-value)", x = 0.92, y = 0.825, gp = gpar(fontsize = 9))
dev.off()


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
