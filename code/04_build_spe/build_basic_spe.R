library("SpatialExperiment")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "V10A27004_A1_Br3874",
        "V10A27004_D1_Br3880",
        "V10A27106_A1_Br3874",
        "V10A27106_B1_Br3854",
        "V10A27106_C1_Br3873",
        "V10A27106_D1_Br3880",
        "V10T31036_A1_Br3874",
        "V10T31036_B1_Br3854",
        "V10T31036_C1_Br3873",
        "V10T31036_D1_Br3880"
    )
)
sample_info$subject <- gsub(".*_", "", sample_info$sample_id)
sample_info$sample_path <- file.path(
    here::here("processed-data", "spaceranger"),
    sample_info$sample_id,
    "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/Visium_IF_AD/blob/master/raw-data/Visium_IF_AD_ITG_MasterExcelSummarySheet.xlsx
donor_info <- data.frame(
    subject = c("Br3854", "Br3873", "Br3880", "Br3874"),
    age = c(65.75, 88.78, 90.47, 73.05),
    sex = c("F", "F", "M", "M"),
    race = "EA/CAUC",
    pmi = c(31.5, 29, 35, 13.5),
    diagnosis = c("AD", "AD", "AD", "Control"),
    rin = c(7, 7.2, 7.1, 7.2),
    BCrating = c("Def AD", "Def AD", "Prob AD", "No AP"),
    braak = c("B3", "B3", "B3", "B2"),
    cerad = c("C3", "C3", "C3", "C0")
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)


## Build basic SPE
Sys.time()
spe <- read10xVisium(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = "lowres",
  load = TRUE
)
Sys.time()
# [1] "2021-10-12 14:49:17 EDT"
# [1] "2021-10-12 14:51:27 EDT"

## Check paths to the targeted sequencing data
stopifnot(all(file.exists(gsub("spaceranger", "spaceranger_targeted", sample_info$sample_path))))

Sys.time()
spe_targeted <- read10xVisium(
  gsub("spaceranger", "spaceranger_targeted", sample_info$sample_path),
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = "lowres",
  load = TRUE
)
Sys.time()
# [1] "2021-10-12 14:51:43 EDT"
# [1] "2021-10-12 14:53:09 EDT"

## This is the case since we didn't use the --target-panel option when
## running spaceranger as described at 
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count
stopifnot(identical(rowData(spe), rowData(spe_targeted)))
stopifnot(identical(colData(spe), colData(spe_targeted)))
spe_targeted$key <- spe$key <- paste0(colnames(spe), '_', spe$sample_id)

## Add the study design info
new_col <- merge(colData(spe), sample_info)
## Fix order
new_col <- new_col[match(spe$key, new_col$key), ]
stopifnot(identical(new_col$key, spe$key))
rownames(new_col) <- rownames(colData(spe))
colData(spe) <- colData(spe_targeted) <- new_col[, -which(colnames(new_col) == "sample_path")]

## Add some information used by spatialLIBD
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
spe_targeted$sum_umi <- colSums(counts(spe_targeted))
spe_targeted$sum_gene <- colSums(counts(spe_targeted) > 0)

## Read in the gene information from the annotation GTF file
gtf <-
  rtracklayer::import(
    "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(spe), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_id, "; ", rowData(spe)$gene_name)
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

## Again for targeted
rowRanges(spe_targeted) <- rowRanges(spe)
spe_targeted$expr_chrM <- colSums(counts(spe_targeted)[is_mito, , drop = FALSE])
spe_targeted$expr_chrM_ratio <- spe_targeted$expr_chrM / spe_targeted$sum_umi


## Read in cell counts and segmentation results
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
    file <- here("processed-data", "spaceranger", sampleid, "outs", "spatial", "tissue_spot_counts.csv")
    if(!file.exists(file)) return(NULL)
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    x
})
## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <- segmentations[segmentation_match, - which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key"))]
colData(spe) <- cbind(colData(spe), segmentation_info)
colData(spe_targeted) <- cbind(colData(spe_targeted), segmentation_info)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8748
length(no_expr) / nrow(spe) * 100
# [1] 23.90099
spe <- spe[-no_expr, ]

no_expr <- which(rowSums(counts(spe_targeted)) == 0)
length(no_expr)
# [1] 13116
length(no_expr) / nrow(spe_targeted) * 100
# [1] 35.83509
spe_targeted <- spe_targeted[-no_expr, ]

## For visualizing this later with spatialLIBD
stopifnot(identical(spatialData(spe), spatialData(spe_targeted)))
spe_targeted$overlaps_tissue <- spe$overlaps_tissue <- factor(ifelse(spatialData(spe)$in_tissue, "in", "out"))

## Size in Mb
lobstr::obj_size(spe) / 1024^2
# 659.3923
lobstr::obj_size(spe_targeted) / 1024^2
# 233.0216

## Save with and without dropping spots outside of the tissue
spe_raw <- spe
spe_raw_targeted <- spe_targeted

dir.create(here::here("processed-data", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "spe", "spe_raw.Rdata"))
save(spe_raw_targeted, file = here::here("processed-data", "spe", "spe_raw_targeted.Rdata"))

## Now drop the spots outside the tissue
spe <- spe_raw[, spatialData(spe_raw)$in_tissue]
dim(spe)
# [1] 27853 38287
## Remove spots without counts
# spe <- spe[, -which(colSums(counts(spe)) == 0)]
# dim(spe)

spe_targeted <- spe_raw_targeted[, spatialData(spe_raw_targeted)$in_tissue]
dim(spe_targeted)
# [1] 23485 38287
## Remove spots without counts
# spe_targeted <- spe_targeted[, -which(colSums(counts(spe_targeted)) == 0)]
# dim(spe_targeted)

lobstr::obj_size(spe) / 1024^2
# 540.7927
lobstr::obj_size(spe_targeted) / 1024^2
# 197.1583

save(spe, file = here::here("processed-data", "spe", "spe.Rdata"))
save(spe_targeted, file = here::here("processed-data", "spe", "spe_targeted.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.1 Patched (2021-08-13 r80752)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-10-12
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  beachmat               2.9.1    2021-08-11 [2] Bioconductor
#  Biobase              * 2.53.0   2021-05-19 [2] Bioconductor
#  BiocGenerics         * 0.39.2   2021-08-18 [1] Bioconductor
#  BiocIO                 1.3.0    2021-05-19 [2] Bioconductor
#  BiocParallel           1.27.17  2021-10-10 [1] Bioconductor
#  Biostrings             2.61.2   2021-08-04 [2] Bioconductor
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  cli                    3.0.1    2021-07-17 [2] CRAN (R 4.1.0)
#  colorout               1.2-2    2021-09-13 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  DelayedArray           0.19.4   2021-09-23 [1] Bioconductor
#  DelayedMatrixStats     1.15.2   2021-08-05 [2] Bioconductor
#  digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.1)
#  DropletUtils           1.13.4   2021-09-19 [1] Bioconductor
#  edgeR                  3.35.0   2021-05-19 [2] Bioconductor
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb         * 1.29.8   2021-09-05 [1] Bioconductor
#  GenomeInfoDbData       1.2.6    2021-05-21 [2] Bioconductor
#  GenomicAlignments      1.29.0   2021-05-19 [2] Bioconductor
#  GenomicRanges        * 1.45.0   2021-05-19 [2] Bioconductor
#  ggplot2                3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array              1.21.0   2021-05-19 [2] Bioconductor
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.1)
#  htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                 1.6.1    2021-05-07 [2] CRAN (R 4.1.0)
#  IRanges              * 2.27.2   2021-08-18 [1] Bioconductor
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  later                  1.2.0    2021-04-23 [2] CRAN (R 4.1.0)
#  lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.1)
#  lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  limma                  3.49.4   2021-08-08 [2] Bioconductor
#  lobstr                 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                 2.7.2    2021-05-02 [2] CRAN (R 4.1.0)
#  magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.1)
#  MatrixGenerics       * 1.5.4    2021-08-26 [1] Bioconductor
#  matrixStats          * 0.61.0   2021-09-17 [1] CRAN (R 4.1.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                 1.6.2    2021-07-29 [2] CRAN (R 4.1.0)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                2.10.1   2020-08-26 [2] CRAN (R 4.1.0)
#  R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
#  Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                  1.98-1.5 2021-09-17 [1] CRAN (R 4.1.1)
#  restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                  2.37.0   2021-05-19 [2] Bioconductor
#  rhdf5filters           1.5.0    2021-05-19 [2] Bioconductor
#  Rhdf5lib               1.15.2   2021-07-01 [2] Bioconductor
#  rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                  0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
#  rmote                  0.3.4    2021-09-13 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools              2.9.1    2021-06-17 [2] Bioconductor
#  rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rtracklayer          * 1.53.1   2021-08-13 [1] Bioconductor
#  S4Vectors            * 0.31.5   2021-10-01 [1] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scuttle                1.3.1    2021-08-05 [1] Bioconductor
#  servr                  0.23     2021-08-11 [1] CRAN (R 4.1.1)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.1.1)
#  SingleCellExperiment * 1.15.1   2021-05-21 [2] Bioconductor
#  sparseMatrixStats      1.5.2    2021-08-05 [2] Bioconductor
#  SpatialExperiment    * 1.3.4    2021-08-24 [1] Bioconductor
#  SummarizedExperiment * 1.23.5   2021-10-05 [1] Bioconductor
#  tibble                 3.1.4    2021-08-25 [1] CRAN (R 4.1.1)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  xfun                   0.25     2021-08-06 [2] CRAN (R 4.1.1)
#  XML                    3.99-0.8 2021-09-17 [1] CRAN (R 4.1.1)
#  XVector                0.33.0   2021-05-19 [2] Bioconductor
#  yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc               1.39.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lcollado/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
