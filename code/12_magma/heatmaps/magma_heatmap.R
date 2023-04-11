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

magmaStats <- list()

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(k)
gene_set <- c("50", "100", "200", "fdr")
print(paste0("Geneset: ", gene_set[k]))


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


#### reshape to long ####
magmaStats_long <- reshape2::melt(magmaStats_wide[, 1:5])
colnames(magmaStats_long)[3:4] <- c("GWAS", "P")
dim(magmaStats_long)

## Calculate Global FDR and Bonf, i.e. across all Pathologies X Disorders ##
magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
magmaStats_long$P.adj.bonf <- p.adjust(magmaStats_long$P, "bonf")


table(p.adjust(magmaStats_long$P, "fdr") < 0.1)

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



## magmaStats_long.beta##
magmaStats_long.beta <- reshape2::melt(magmaStats_wide.beta)
colnames(magmaStats_long.beta)[3:4] <- c("GWAS", "Beta")
magmaStats_long.beta

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

print("magmaStats_wide")
magmaStats_wide

print("magmaStats_long")
magmaStats_long

## Print to CSV ===
write.csv(magmaStats_long,
    file = here(
        "code", "12_magma", "heatmaps",
        paste0("Global_FDRs_Bonf_3GWAS_x_7_Pathologies_", gene_set[k], ".csv")
    ),
    row.names = F, quote = F
)

write.csv(magmaStats_wide,
    file = here(
        "code", "12_magma", "heatmaps",
        paste0("GWAS_wise_FDRs_Bonf_3GWAS_x_7_Pathologies_", gene_set[k], ".csv")
    ),
    row.names = F, quote = F
)

head(magmaStats_long)

#### Make heatmap ####
midpoint <- function(x) x[-length(x)] + diff(x) / 2

MAGMAplot <- function(region, Pthresh, fdrThresh, ...) {
    #### Set up -log10(p_value) ####
    wide_p <- sapply(magmaStats[[region]], function(x) {
        cbind(-log10(x$P))
    })
    rownames(wide_p) <- magmaStats[[region]][[1]]$VARIABLE
    wide_p[wide_p > Pthresh] <- Pthresh
    wide_p <- round(wide_p[rev(sort(rownames(wide_p))), ], 3)


    #### Set up betas ####
    wide_beta <- sapply(magmaStats[[region]], function(x) {
        cbind(x$BETA)
    })
    rownames(wide_beta) <- magmaStats[[region]][[1]]$VARIABLE
    wide_beta <- round(wide_beta[rev(sort(rownames(wide_beta))), ], 2)

    #### Use empirical cutoff (independent p=0.05) for printing betas ####
    wide_beta[wide_p < -log10(0.05)] <- ""

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
        col = mypal, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...
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
dir_output <- dir.create(here("plots", "12_magma"))
pdf(here("plots", "12_magma", paste0("heatmap_", "top_", gene_set[k], "_genes", ".pdf")), w = 6)
par(mar = c(8.5, 7.5, 6, 1), cex.axis = 1.0, cex.lab = 0.5)
MAGMAplot(region = "ITC", Pthresh = 12, fdrThresh = 1) ## set to 1 because didn't want to make a cutoff based on FDR/Bonf
abline(v = 7, lwd = 3)
text(x = 1.5, y = 185, "MAGMA gene set analyses: ITC pathology types", xpd = TRUE, cex = 1.5, font = 2)
text(x = 1.5, y = 175, "(Betas for empirically significant associations, p < 0.05)", xpd = TRUE, cex = 1, font = 1)
grid::grid.text(label = "-log10(p-value)", x = 0.92, y = 0.825, gp = gpar(fontsize = 9))
dev.off()


## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
