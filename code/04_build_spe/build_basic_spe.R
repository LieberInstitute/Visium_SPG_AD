library("SpatialExperiment")
library("spatialLIBD")
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
spe <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE
)
Sys.time()
# [1] "2021-11-10 17:57:41 EST"
# 2021-11-10 17:57:42 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2021-11-10 18:00:33 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2021-11-10 18:00:37 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2021-11-10 18:00:40 rtracklayer::import: reading the reference GTF file
# 2021-11-10 18:02:15 adding gene information to the SPE object
# 2021-11-10 18:02:15 adding information used by spatialLIBD
# [1] "2021-11-10 18:02:22 EST"

## Check paths to the targeted sequencing data
stopifnot(all(file.exists(gsub("spaceranger", "spaceranger_targeted", sample_info$sample_path))))

Sys.time()
spe_targeted <- read10xVisiumWrapper(
  gsub("spaceranger", "spaceranger_targeted", sample_info$sample_path),
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres"),
  load = TRUE
)
Sys.time()
# [1] "2021-10-12 14:51:43 EDT"
# [1] "2021-10-12 14:53:09 EDT"


## Add images created by Madhavi Tippani
spe <- add_images(
    spe = spe,
    image_dir = here("processed-data", "Images", "spatialLIBD_images"),
    image_pattern = "Abeta_lowres",
    image_id_current = "lowres"
)
spe <- add_images(
    spe = spe,
    image_dir = here("processed-data", "Images", "spatialLIBD_images"),
    image_pattern = "Abeta_hires",
    image_id_current = "hires"
)

## Update the imData in the targeted sequencing
imgData(spe_targeted) <- imgData(spe)

## This is the case since we didn't use the --target-panel option when
## running spaceranger as described at
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count
stopifnot(identical(rowData(spe), rowData(spe_targeted)))
stopifnot(identical(colData(spe), colData(spe_targeted)))

## Add the study design info
new_col <- merge(colData(spe), sample_info)
## Fix order
new_col <- new_col[match(spe$key, new_col$key), ]
stopifnot(identical(new_col$key, spe$key))
rownames(new_col) <- rownames(colData(spe))
colData(spe_targeted) <- colData(spe) <- new_col[, -which(colnames(new_col) == "sample_path")]

## Read in cell counts and segmentation results
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
    file <- here("processed-data", "spaceranger", sampleid, "outs", "spatial", "tissue_spot_counts.csv")
    if(!file.exists(file)) return(NULL)
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
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

## Save with and without dropping spots outside of the tissue
spe_raw <- spe
spe_raw_targeted <- spe_targeted

dir.create(here::here("processed-data", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "spe", "spe_raw.Rdata"))
save(spe_raw_targeted, file = here::here("processed-data", "spe", "spe_raw_targeted.Rdata"))

## Size in Mb
lobstr::obj_size(spe_raw) / 1024^2
# 1,595.841
lobstr::obj_size(spe_raw_targeted) / 1024^2
# 233.0216

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

# ─ Session info  ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  hash: camera, steaming bowl, prince
#
#  setting  value
#  version  R version 4.1.2 Patched (2021-11-04 r81138)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-11-10
#  pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date (UTC) lib source
#  AnnotationDbi            1.56.1   2021-10-29 [2] Bioconductor
#  AnnotationHub            3.2.0    2021-10-26 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.2)
#  beachmat                 2.10.0   2021-10-26 [2] Bioconductor
#  beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.2)
#  benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.2)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.2)
#  Biobase                * 2.54.0   2021-10-26 [2] Bioconductor
#  BiocFileCache            2.2.0    2021-10-26 [2] Bioconductor
#  BiocGenerics           * 0.40.0   2021-10-26 [2] Bioconductor
#  BiocIO                   1.4.0    2021-10-26 [2] Bioconductor
#  BiocManager              1.30.16  2021-06-15 [2] CRAN (R 4.1.2)
#  BiocNeighbors            1.12.0   2021-10-26 [1] Bioconductor
#  BiocParallel             1.28.0   2021-10-26 [2] Bioconductor
#  BiocSingular             1.10.0   2021-10-26 [1] Bioconductor
#  BiocVersion              3.14.0   2021-05-19 [2] Bioconductor
#  Biostrings               2.62.0   2021-10-26 [2] Bioconductor
#  bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  blob                     1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
#  bslib                    0.3.1    2021-10-06 [2] CRAN (R 4.1.2)
#  cachem                   1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
#  callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
#  cli                      3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
#  codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
#  colorout                 1.2-2    2021-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.2)
#  cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.2)
#  crayon                   1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
#  curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
#  data.table               1.14.2   2021-09-27 [2] CRAN (R 4.1.2)
#  DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
#  DelayedArray             0.20.0   2021-10-26 [2] Bioconductor
#  DelayedMatrixStats       1.16.0   2021-10-26 [2] Bioconductor
#  desc                     1.4.0    2021-09-28 [2] CRAN (R 4.1.2)
#  digest                   0.6.28   2021-09-23 [2] CRAN (R 4.1.2)
#  dockerfiler              0.1.4    2021-09-03 [1] CRAN (R 4.1.2)
#  doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)
#  dplyr                    1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.2)
#  DropletUtils             1.14.1   2021-11-08 [1] Bioconductor
#  DT                       0.19     2021-09-02 [2] CRAN (R 4.1.2)
#  edgeR                    3.36.0   2021-10-26 [2] Bioconductor
#  ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  ExperimentHub            2.2.0    2021-10-26 [2] Bioconductor
#  fansi                    0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   13.3     2021-10-30 [2] CRAN (R 4.1.2)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics                 0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
#  GenomeInfoDb           * 1.30.0   2021-10-26 [2] Bioconductor
#  GenomeInfoDbData         1.2.7    2021-11-01 [2] Bioconductor
#  GenomicAlignments        1.30.0   2021-10-26 [2] Bioconductor
#  GenomicRanges          * 1.46.0   2021-10-26 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.2)
#  ggplot2                  3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  glue                     1.5.0    2021-11-07 [2] CRAN (R 4.1.2)
#  golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.2)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.22.0   2021-10-26 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.2)
#  htmltools                0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
#  htmlwidgets              1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
#  httpuv                   1.6.3    2021-09-09 [2] CRAN (R 4.1.2)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  interactiveDisplayBase   1.32.0   2021-10-26 [2] Bioconductor
#  IRanges                * 2.28.0   2021-10-26 [2] Bioconductor
#  irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  KEGGREST                 1.34.0   2021-10-26 [2] Bioconductor
#  knitr                    1.36     2021-09-29 [2] CRAN (R 4.1.2)
#  later                    1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
#  lattice                  0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
#  limma                    3.50.0   2021-10-26 [2] Bioconductor
#  lobstr                 * 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
#  magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  maps                     3.4.0    2021-09-25 [2] CRAN (R 4.1.2)
#  Matrix                   1.3-4    2021-06-01 [3] CRAN (R 4.1.2)
#  MatrixGenerics         * 1.6.0    2021-10-26 [2] Bioconductor
#  matrixStats            * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
#  memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
#  mime                     0.12     2021-09-28 [2] CRAN (R 4.1.2)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
#  pkgbuild                 1.2.0    2020-12-15 [2] CRAN (R 4.1.0)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.2.3    2021-10-13 [2] CRAN (R 4.1.2)
#  plotly                   4.10.0   2021-10-09 [2] CRAN (R 4.1.2)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome               1.3.1    2021-07-16 [1] CRAN (R 4.1.2)
#  prettyunits              1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
#  processx                 3.5.2    2021-04-30 [2] CRAN (R 4.1.0)
#  promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  ps                       1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
#  purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                  2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
#  R6                       2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
#  rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                     1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                    1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
#  remotes                  2.4.1    2021-09-29 [2] CRAN (R 4.1.2)
#  restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                    2.38.0   2021-10-26 [2] Bioconductor
#  rhdf5filters             1.6.0    2021-10-26 [2] Bioconductor
#  Rhdf5lib                 1.16.0   2021-10-26 [2] Bioconductor
#  rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                    0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
#  rmote                    0.3.4    2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.1.2    2021-09-08 [2] CRAN (R 4.1.2)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools                2.10.0   2021-10-26 [2] Bioconductor
#  RSQLite                  2.2.8    2021-08-21 [2] CRAN (R 4.1.2)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.2)
#  rtracklayer            * 1.54.0   2021-10-26 [2] Bioconductor
#  S4Vectors              * 0.32.2   2021-11-07 [2] Bioconductor
#  sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
#  ScaledMatrix             1.2.0    2021-10-26 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                   1.22.0   2021-10-26 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.2)
#  scuttle                  1.4.0    2021-10-26 [1] Bioconductor
#  servr                    0.23     2021-08-11 [1] CRAN (R 4.1.2)
#  sessioninfo            * 1.2.1    2021-11-02 [2] CRAN (R 4.1.2)
#  shiny                    1.7.1    2021-10-02 [2] CRAN (R 4.1.2)
#  shinyWidgets             0.6.2    2021-09-17 [1] CRAN (R 4.1.2)
#  SingleCellExperiment   * 1.16.0   2021-10-26 [2] Bioconductor
#  spam                     2.7-0    2021-06-25 [2] CRAN (R 4.1.0)
#  sparseMatrixStats        1.6.0    2021-10-26 [2] Bioconductor
#  SpatialExperiment      * 1.4.0    2021-10-26 [1] Bioconductor
#  spatialLIBD            * 1.7.3    2021-11-10 [1] Github (LieberInstitute/spatialLIBD@771d2f7)
#  stringi                  1.7.5    2021-10-04 [2] CRAN (R 4.1.2)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.24.0   2021-10-26 [2] Bioconductor
#  testthat                 3.1.0    2021-10-04 [2] CRAN (R 4.1.2)
#  tibble                   3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
#  tidyr                    1.1.4    2021-09-27 [2] CRAN (R 4.1.2)
#  tidyselect               1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  usethis                  2.1.3    2021-10-27 [2] CRAN (R 4.1.2)
#  utf8                     1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.2)
#  viridis                  0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
#  viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  withr                    2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  xfun                     0.28     2021-11-04 [2] CRAN (R 4.1.2)
#  XML                      3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
#  xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
#  xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                  0.34.0   2021-10-26 [2] Bioconductor
#  yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc                 1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
