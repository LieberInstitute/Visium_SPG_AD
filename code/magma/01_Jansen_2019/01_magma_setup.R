#borrowed from https://github.com/lmweber/locus-c/blob/main/code/magma/MAGMA_setup.R


#using R 4.1.x since jaffelab requires a stable version
library(sgejobs)
library(rtracklayer)
library(GenomicRanges)
#BiocManager::install("liftOver")
library(liftOver)
library(jaffelab)
library(sessioninfo)
library(here)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
here()
#[1] "/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD"



####Create gene location map for SNPs ####

## get GTF
gtf = import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

#[1] 2765969
gtf = gtf[gtf$type == "gene"]
length(gtf)
#36601

#### what does a GTF row look like? ####
# > gtf[1,]
# GRanges object with 1 range and 24 metadata columns:
#     seqnames      ranges strand |   source     type     score     phase
# <Rle>   <IRanges>  <Rle> | <factor> <factor> <numeric> <integer>
#     [1]     chr1 29554-31109      + |   HAVANA     gene        NA      <NA>
#     gene_id gene_version   gene_type   gene_name       level
# <character>  <character> <character> <character> <character>
#     [1] ENSG00000243485            5      lncRNA MIR1302-2HG           2
# hgnc_id         tag          havana_gene transcript_id
# <character> <character>          <character>   <character>
#     [1]  HGNC:52482  ncRNA_host OTTHUMG00000000959.2          <NA>
#     transcript_version transcript_type transcript_name
# <character>     <character>     <character>
#     [1]               <NA>            <NA>            <NA>
#     transcript_support_level havana_transcript exon_number     exon_id
# <character>       <character> <character> <character>
#     [1]                     <NA>              <NA>        <NA>        <NA>
#     exon_version  protein_id      ccdsid         ont
# <character> <character> <character> <character>
#     [1]         <NA>        <NA>        <NA>        <NA>
#     -------
#     seqinfo: 40 sequences from an unspecified genome; no seqlengths



#extract relevant columns from GTF
annoTab.full <- as.data.frame(ranges(gtf))
annoTab.full$names <- mcols(gtf)$gene_id
annoTab.full$symbol <- mcols(gtf)$gene_name
annoTab.full$seqlevel <- as.character(seqnames(gtf))
annoTab.full$strand <- as.character(strand(gtf))
annoTab.full <- annoTab.full[ ,c("names", "seqlevel", "start", "end", "strand",
                                 "symbol")]
table(annoTab.full$seqlevel)

#### ####

# chr1      chr10      chr11      chr12      chr13      chr14      chr15
# 3410       1394       2065       1928        788       1474       1267
# chr16      chr17      chr18      chr19       chr2      chr20      chr21
# 1649       1992        781       2027       2541        964        555
# chr22       chr3       chr4       chr5       chr6       chr7       chr8
# 901       1893       1533       1811       1827       1686       1495
# chr9       chrM       chrX       chrY GL000009.2 GL000194.1 GL000195.1
# 1319         13       1148        111          1          2          2
# GL000205.2 GL000213.1 GL000218.1 GL000219.1 KI270711.1 KI270713.1 KI270721.1
# 1          1          1          1          1          2          1
# KI270726.1 KI270727.1 KI270728.1 KI270731.1 KI270734.1
# 2          4          6          1          3
dir.create(here("code","magma","annotation"))

#### ####

write.table(annoTab.full, file=here("code","magma","annotation",
                                    "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"),
            sep="\t",
            row.names=F, col.names=F, quote=F)


## Define all expressed genes

marker_list <- read.table(here("code","magma", "01_Jansen_2019","pvalues_top_100.txt"),
                      header = TRUE)



### Lift coordinates from GRCh38 (hg38) > GRCh37 (hg19) ===
# (because GWAS are still run w/ summary stats in hg19 coords)
# - Also only use those for expressed genes - otherwise MTC penalizes stats
gene_df = read.delim(here("code","magma","annotation",
                          "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"),
                     header=FALSE)
colnames(gene_df)= c("GeneID", "Chr", "Start", "End", "Strand", "Symbol")

# add 'chr' prefix

gr = makeGRangesFromDataFrame(gene_df, keep=TRUE)
names(gr) = gr$GeneID

## Lift to hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
lifted_list = range(liftOver(gr, ch))
# Discarding unchained sequences: GL000009.2, GL000194.1, GL000195.1,
# GL000205.2, GL000213.1, GL000218.1, GL000219.1, KI270711.1, KI270713.1,
# KI270721.1, KI270726.1, KI270727.1, KI270728.1, KI270731.1, KI270734.1


table(lengths(lifted_list))
# 0     1     2     3     4     5     6     7    11    14    16    17    21  31
# 78 36377   116    12     6     3     1     2     1     1     1     1     1 1

lifted_list = lifted_list[lengths(lifted_list) == 1]


lifted = unlist(lifted_list)
lifted_df = as.data.frame(lifted)
lifted_df$Symbol = gene_df$Symbol[match(rownames(lifted_df), gene_df$GeneID)]

# Check:
lifted_df[lifted_df$Symbol=="GABRQ", ]
#                seqnames     start       end width strand Symbol
#ENSG00000268089     chrX 151806349 151826007 19659      +  GABRQ
#   ^ good - in hg19, has a different EnsemblID ("ENSG00000147402")

# And in hg38 coords:
gene_df[gene_df$Symbol=="GABRQ", ]
#               GeneID  Chr     Start       End Strand Symbol
#36340 ENSG00000268089 chrX 152637895 152657542      +  GABRQ


## Rearrange cols of lifted_df / match format ===
table(marker_list$Gene %in% rownames(lifted_df))
# FALSE  TRUE
# 1   699

head(lifted_df, n=3)
#                 seqnames  start    end width strand      Symbol
# ENSG00000243485     chr1  29554  31109  1556      + MIR1302-2HG
# ENSG00000237613     chr1  34554  36081  1528      -     FAM138A
# ENSG00000186092     chr1  65419  71585  6167      +       OR4F5

head(gene_df, n=3)
#            GeneID  Chr  Start    End Strand      Symbol
# 1 ENSG00000243485 chr1  29554  31109      + MIR1302-2HG
# 2 ENSG00000237613 chr1  34554  36081      -     FAM138A
# 3 ENSG00000186092 chr1  65419  71585      +       OR4F5     <- match this format for MAGMA input

lifted_df$width <- NULL # don't need
lifted_df$GeneID <- rownames(lifted_df)
lifted_df <- lifted_df[ ,c(6,1:5)]

# Remove the 'chr' prefix
lifted_df$seqnames <- ss(as.character(lifted_df$seqnames),"chr",2)

# Finally subset for those expressed genes
lifted_df <- lifted_df[intersect(rownames(lifted_df), marker_list$Gene), ]

## Write out for MAGMA SNP annotation step
write.table(lifted_df, file=here("code","magma","annotation",
                                 "GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc"), sep="\t",
            row.names=F, col.names=F, quote=F)




### GWAS SNP location files ====================================================
# ** As per manual, the column order should be: SNP ID, chromosome, bp position

## GWAS for AD (PGC-ALZ, IGAP, ADSP, UKB meta-meta-analysis: Jansen, et al.2019) ===


sumStats.AD <- read.table("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt", header=T)
class(sumStats.AD)
dim(sumStats.AD)
#[1] 13367299       14
head(sumStats.AD)
unique(sumStats.AD$CHR)  # 1:22 (no MT genes)
snploc.AD <- sumStats.AD[ ,c("SNP", "CHR", "BP")]
write.table(snploc.AD, file=here("code","magma","GWAS_Results","Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc"),
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".AD"))



#




## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# [1] "2022-05-29 16:56:28 EDT"
proc.time()
#    user   system  elapsed
#  46.493    1.946 1202.705
options(width = 120)
session_info()
