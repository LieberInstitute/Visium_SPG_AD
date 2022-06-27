# GWAS for FTD () ===

# library(sgejobs)
# sgejobs::job_single(
#     "01_magma_setup_FTD",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     command = "Rscript 01_percent_reads_assigned.R",
#     create_logdir = TRUE
# )


library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(biomaRt)
#BiocManager::install("liftOver")
library(liftOver)
library(jaffelab)
library(sessioninfo)
library(here)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
here()

#sumStats.FTD <- read.table(here("raw-data","magma_GWAS_files","FTD_GWAS_META.txt"), header=T)
sumStats.FTD <- read.table(here("raw-data","magma_GWAS_files","FTD_GWAS_META.txt"), header=T)
class(sumStats.FTD)
dim(sumStats.FTD)
#[1]  6026384       8
head(sumStats.FTD)
# marker Allele1 Allele2   beta1     SE  pValue chr        Bp
# chr10:100004799       a       c  0.3528 0.2808 0.20890  10 100004799
##Similar to Matt's PD GWAS, this only has the SNP location information but not the
##rsIDs.
unique(sumStats.FTD$chr)  # 1:22
# [1] 10  1 11 12 13 14 15 16 17 18 19 20  2 21 22  3  4  5  6  7  8  9

###library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snpcount(snps)

snps_38 <- SNPlocs.Hsapiens.dbSNP151.GRCh38
snpcount(snps_38)
# Try it with the smallest autosome:
chr21_snps_38 <- snpsBySeqname(snps_38, "21")
class(chr21_snps_38)
head(chr21_snps_38)

chr21_snps_38 <- as.data.frame(chr21_snps_38)
chr21_snps_38$chr.bp <- paste0("chr",chr21_snps_38$seqnames,":",chr21_snps_38$pos)

table(sumStats.FTD$chr == "21")

sumStats.FTD.chr21 <- sumStats.FTD[sumStats.FTD$chr=="21", ]


table(sumStats.FTD.chr21$marker %in% chr21_snps_38$chr.bp)
# FALSE  TRUE
# 70014 15597

sumStats.FTD.keep <- data.frame()

temp.df <- sumStats.FTD.chr21 |> dplyr::filter(marker %in%  chr21_snps$chr.bp)
temp.df$rsID <- chr21_snps$RefSNP_id[match(temp.df$marker, chr21_snps$chr.bp)]

# Rbind
sumStats.FTD.keep <- rbind(sumStats.FTD.keep, temp.df)

for(i in seqnames(snps)){
    cat("Querying chr: ",i,"...\n")
    temp.snps <- as.data.frame(snpsBySeqname(snps, i))
    temp.snps$chr.bp <- paste0("chr",temp.snps$seqnames,":",temp.snps$pos)

    # Subset sumStats to quantify % intersecting per chromosome
    sumStats.temp <- sumStats.FTD[sumStats.FTD$chr==i, ]

    cat(paste0("\tPercent of summary statistics SNPs in chr:",i," with rsIDs:\n"))
    print(table(sumStats.temp$marker %in% temp.snps$chr.bp)["TRUE"] / nrow(sumStats.temp) * 100)
    cat("\n")

    temp.df <- sumStats.temp[intersect(sumStats.temp$marker, temp.snps$chr.bp), ]
    temp.df$rsID <- temp.snps$RefSNP_id[match(temp.df$marker, temp.snps$chr.bp)]

    # Rbind
    sumStats.FTD.keep <- rbind(sumStats.FTD.keep, temp.df)
}


# Total SNPs in summary stats with rsIDs:
dim(sumStats.FTD.keep)
(nrow(sumStats.FTD.keep) / nrow(sumStats.FTD)) * 100

snploc.FTD <- sumStats.FTD.keep[ ,c("rsID", "chr", "Bp")]
write.table(snploc.FTD, file=here("code","magma","02_lancet_2014","FTD_Lancet2014.snploc"),
            sep="\t", col.names=T, row.names=F, quote=F)

## Create an 'Neff' using METAL's recommended computation for meta-GWAS (https://doi.org/10.1093/bioinformatics/btq340)
#      instead of sum(N_cases, N_controls)
#sumStats.FTD.keep$N_effective <- 4/(1/sumStats.PD.keep$N_cases + 1/sumStats.PD.keep$N_controls)

# Save
write.table(sumStats.FTD.keep, file=here("code","magma","02_lancet_2014",
                                        "FTD-IFGC-and-rsID-ADDED.tab"),
            sep="\t", col.names=T, row.names=F, quote=F)

rm(list=ls(pattern=".PD"))
