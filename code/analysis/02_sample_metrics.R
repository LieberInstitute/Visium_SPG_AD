## Automatically style the code in this script:
styler::style_file(here("code", "analysis", "02_sample_metrics.R"),
    transformers = biocthis::bioc_style()
)


library("ggplot2")
library("here")
library("ggpubr")
library("sessioninfo")


sample_names <-
    c(
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
dir_outputs <- c(
    here("raw-data", "10x_files", "Lieber_Transfer"),
    here("raw-data", "10x_files", "Lieber_Transfer_10x_Alignments")
)

metrics_csvs <- unlist(lapply(dir_outputs, function(x) {
    file.path(x, sample_names)
}))
stopifnot(all(file.exists(metrics_csvs)))

df_metrics_all <- NULL
for (i in metrics_csvs) {
    dir_csv <- file.path(i, "metrics_summary_csv.csv")
    df_metrics <- read.csv(dir_csv, header = TRUE)
    df_metrics_all <- rbind(df_metrics_all, df_metrics)
}

# save
sample_metrics <- df_metrics_all

## Overwrite the 10x Sample IDs by ours
sample_metrics$Sample.ID <- basename(metrics_csvs)
sample_metrics$Alignment <- rep(c("Abby", "10x"), each = 10)

dir.create(here("processed-data", "10x_checks"), showWarnings = FALSE)
save(sample_metrics,
    file = here("processed-data", "10x_checks", "sample_metrics.Rdata")
)
write.csv(sample_metrics,
    file = here("processed-data", "10x_checks", "sample_metrics.csv")
)


dir.create(here("plots", "10x_checks"), showWarnings = FALSE)
pdf(
    here("plots", "10x_checks", "spaceranger_metrics_by_number_of_reads.pdf"),
    useDingbats = FALSE,
    width = 10
)
ggplot(
    sample_metrics,
    aes(x = Mean.Reads.per.Spot, y = Number.of.Reads / 1e6)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20) +
    facet_grid(~
    Alignment)

ggplot(
    sample_metrics,
    aes(x = Median.Genes.per.Spot, y = Number.of.Reads / 1e6, color = Alignment)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)
ggplot(
    sample_metrics,
    aes(x = Total.Genes.Detected, y = Number.of.Reads / 1e6, color = Alignment)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)
ggplot(
    sample_metrics,
    aes(x = Median.UMI.Counts.per.Spot, y = Number.of.Reads / 1e6, color = Alignment)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)

dev.off()



## Compare fiducial alignments: made by Abby or by 10x Genomics
## Used https://rpkgs.datanovia.com/ggpubr/reference/ggpaired.html
pdf(
    here("plots", "10x_checks", "alignment_check_abby_vs_10x.pdf"),
    useDingbats = FALSE,
    width = 10
)
for(i in colnames(sample_metrics)[-c(1, ncol(sample_metrics))]) {
    p <- ggpaired(sample_metrics, id = "Sample.ID", x = "Alignment", y = i, xlab = "Fiducial Alignment", ylab = i, fill = "Alignment", palette = "Dark2", line.color = "gray", ggtheme = theme_pubr(base_size = 30))
    print(p)
}
dev.off()

## Match with spatialDLPFC metrics
load(here("raw-data", "spatialDLPFC", "shared_metrics.Rdata"), verbose = TRUE)
shared_metrics$Sample.ID <- rownames(shared_metrics)

## Keep only the alignments by Abby for now
tmp <- subset(sample_metrics, Alignment == "Abby")
tmp$Alignment <- NULL
tmp$study <- "VisiumIF"

## Columns that aren't shared
colnames(tmp)[!colnames(tmp) %in% colnames(shared_metrics)]
# [1] "Number.of.Spots.Under.Tissue"         "Mean.Reads.Under.Tissue.per.Spot"
# [3] "Fraction.of.Spots.Under.Tissue"       "Valid.UMIs"
# [5] "Fraction.Reads.in.Spots.Under.Tissue"
colnames(shared_metrics)[!colnames(shared_metrics) %in% colnames(tmp)]
# character(0)

## Subset to common columns
cols_shared <- colnames(tmp)[colnames(tmp) %in% colnames(shared_metrics)]
tmp <- tmp[, cols_shared]
shared_metrics <- shared_metrics[, cols_shared]

## Put on the same scale
tmp[, sapply(tmp, max) < 1] <- round(tmp[, sapply(tmp, max) < 1] * 100, 1)

## Combine
all_metrics <- rbind(tmp, shared_metrics)
row.names(all_metrics) <- NULL

## Fix study
all_metrics$study[all_metrics$study == "current"] <- "spatialDLPFC"
all_metrics$study[all_metrics$study == "pilot"] <- "spatialLIBD"

## Save for later
save(all_metrics,
    file = here("processed-data", "10x_checks", "all_metrics.Rdata")
)
write.csv(all_metrics,
    file = here("processed-data", "10x_checks", "all_metrics.csv")
)


pdf(
    here("plots", "10x_checks", "cross_study_spaceranger_metrics_by_number_of_reads.pdf"),
    useDingbats = FALSE,
    width = 12
)
ggplot(
    all_metrics,
    aes(x = Mean.Reads.per.Spot, y = Number.of.Reads / 1e6)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20) +
    facet_grid(~
    study)

ggplot(
    all_metrics,
    aes(x = Median.Genes.per.Spot, y = Number.of.Reads / 1e6, color = study)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)
ggplot(
    all_metrics,
    aes(x = Total.Genes.Detected, y = Number.of.Reads / 1e6, color = study)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)
ggplot(
    all_metrics,
    aes(x = Median.UMI.Counts.per.Spot, y = Number.of.Reads / 1e6, color = study)
) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw(base_size = 20)

dev.off()


## Compare metrics one at a time across studies
## Used https://rpkgs.datanovia.com/ggpubr/reference/ggboxplot.html
pdf(
    here("plots", "10x_checks", "cross_study_spaceranger_metrics_boxplots.pdf"),
    useDingbats = FALSE,
    width = 10
)
for(i in colnames(all_metrics)[-c(1, ncol(all_metrics))]) {
    set.seed(20210714)
    p <- ggboxplot(all_metrics, x = "study", y = i, color = "study", palette = "Dark2", add = "jitter", shape = "study", label = "Sample.ID", repel = TRUE, font.label = list(size = 5), ggtheme = theme_pubr(base_size = 30))
    print(p)
}
dev.off()



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.0 (2021-05-18)
#  os       macOS Big Sur 11.4
#  system   x86_64, darwin17.0
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2021-07-14
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package      * version    date       lib source
#  abind          1.4-5      2016-07-21 [1] CRAN (R 4.1.0)
#  assertthat     0.2.1      2019-03-21 [1] CRAN (R 4.1.0)
#  backports      1.2.1      2020-12-09 [1] CRAN (R 4.1.0)
#  biocthis       1.2.0      2021-05-19 [1] Bioconductor
#  broom          0.7.7      2021-06-13 [1] CRAN (R 4.1.0)
#  cachem         1.0.5      2021-05-15 [1] CRAN (R 4.1.0)
#  callr          3.7.0      2021-04-20 [1] CRAN (R 4.1.0)
#  car            3.0-10     2020-09-29 [1] CRAN (R 4.1.0)
#  carData        3.0-4      2020-05-22 [1] CRAN (R 4.1.0)
#  cellranger     1.1.0      2016-07-27 [1] CRAN (R 4.1.0)
#  cli            2.5.0      2021-04-26 [1] CRAN (R 4.1.0)
#  colorout       1.2-2      2020-11-03 [1] Github (jalvesaq/colorout@726d681)
#  colorspace     2.0-1      2021-05-04 [1] CRAN (R 4.1.0)
#  crayon         1.4.1      2021-02-08 [1] CRAN (R 4.1.0)
#  curl           4.3.1      2021-04-30 [1] CRAN (R 4.1.0)
#  data.table     1.14.0     2021-02-21 [1] CRAN (R 4.1.0)
#  DBI            1.1.1      2021-01-15 [1] CRAN (R 4.1.0)
#  desc           1.3.0      2021-03-05 [1] CRAN (R 4.1.0)
#  devtools     * 2.4.2      2021-06-07 [1] CRAN (R 4.1.0)
#  digest         0.6.27     2020-10-24 [1] CRAN (R 4.1.0)
#  dplyr          1.0.7      2021-06-18 [1] CRAN (R 4.1.0)
#  ellipsis       0.3.2      2021-04-29 [1] CRAN (R 4.1.0)
#  fansi          0.5.0      2021-05-25 [1] CRAN (R 4.1.0)
#  farver         2.1.0      2021-02-28 [1] CRAN (R 4.1.0)
#  fastmap        1.1.0      2021-01-25 [1] CRAN (R 4.1.0)
#  forcats        0.5.1      2021-01-27 [1] CRAN (R 4.1.0)
#  foreign        0.8-81     2020-12-22 [1] CRAN (R 4.1.0)
#  fs             1.5.0      2020-07-31 [1] CRAN (R 4.1.0)
#  generics       0.1.0      2020-10-31 [1] CRAN (R 4.1.0)
#  ggplot2      * 3.3.4      2021-06-16 [1] CRAN (R 4.1.0)
#  ggpubr       * 0.4.0      2020-06-27 [1] CRAN (R 4.1.0)
#  ggrepel        0.9.1      2021-01-15 [1] CRAN (R 4.1.0)
#  ggsignif       0.6.2      2021-06-14 [1] CRAN (R 4.1.0)
#  glue           1.4.2      2020-08-27 [1] CRAN (R 4.1.0)
#  gtable         0.3.0      2019-03-25 [1] CRAN (R 4.1.0)
#  haven          2.4.1      2021-04-23 [1] CRAN (R 4.1.0)
#  here         * 1.0.1      2020-12-13 [1] CRAN (R 4.1.0)
#  hms            1.1.0      2021-05-17 [1] CRAN (R 4.1.0)
#  labeling       0.4.2      2020-10-20 [1] CRAN (R 4.1.0)
#  lattice        0.20-44    2021-05-02 [1] CRAN (R 4.1.0)
#  lifecycle      1.0.0      2021-02-15 [1] CRAN (R 4.1.0)
#  lubridate      1.7.10     2021-02-26 [1] CRAN (R 4.1.0)
#  magrittr       2.0.1      2020-11-17 [1] CRAN (R 4.1.0)
#  Matrix         1.3-4      2021-06-01 [1] CRAN (R 4.1.0)
#  memoise        2.0.0      2021-01-26 [1] CRAN (R 4.1.0)
#  mgcv           1.8-36     2021-06-01 [1] CRAN (R 4.1.0)
#  munsell        0.5.0      2018-06-12 [1] CRAN (R 4.1.0)
#  nlme           3.1-152    2021-02-04 [1] CRAN (R 4.1.0)
#  openxlsx       4.2.4      2021-06-16 [1] CRAN (R 4.1.0)
#  pillar         1.6.1      2021-05-16 [1] CRAN (R 4.1.0)
#  pkgbuild       1.2.0      2020-12-15 [1] CRAN (R 4.1.0)
#  pkgconfig      2.0.3      2021-03-31 [1] Github (gaborcsardi/pkgconfig@b81ae03)
#  pkgload        1.2.1      2021-04-06 [1] CRAN (R 4.1.0)
#  prettyunits    1.1.1      2020-01-24 [1] CRAN (R 4.1.0)
#  processx       3.5.2      2021-04-30 [1] CRAN (R 4.1.0)
#  prompt         1.0.1      2021-03-12 [1] CRAN (R 4.1.0)
#  ps             1.6.0      2021-02-28 [1] CRAN (R 4.1.0)
#  purrr          0.3.4      2020-04-17 [1] CRAN (R 4.1.0)
#  R6             2.5.0      2020-10-28 [1] CRAN (R 4.1.0)
#  RColorBrewer   1.1-2      2014-12-07 [1] CRAN (R 4.1.0)
#  Rcpp           1.0.6      2021-01-15 [1] CRAN (R 4.1.0)
#  readxl         1.3.1      2019-03-13 [1] CRAN (R 4.1.0)
#  remotes        2.4.0      2021-06-02 [1] CRAN (R 4.1.0)
#  rio            0.5.27     2021-06-21 [1] CRAN (R 4.1.0)
#  rlang          0.4.11     2021-04-30 [1] CRAN (R 4.1.0)
#  rprojroot      2.0.2      2020-11-15 [1] CRAN (R 4.1.0)
#  rstatix        0.7.0      2021-02-13 [1] CRAN (R 4.1.0)
#  rsthemes       0.2.1.9000 2021-02-12 [1] Github (gadenbuie/rsthemes@521572b)
#  rstudioapi     0.13       2020-11-12 [1] CRAN (R 4.1.0)
#  scales         1.1.1      2020-05-11 [1] CRAN (R 4.1.0)
#  sessioninfo  * 1.1.1      2018-11-05 [1] CRAN (R 4.1.0)
#  stringi        1.6.2      2021-05-17 [1] CRAN (R 4.1.0)
#  styler         1.4.1      2021-03-30 [1] CRAN (R 4.1.0)
#  suncalc        0.5.0      2019-04-03 [1] CRAN (R 4.1.0)
#  testthat     * 3.0.3      2021-06-16 [1] CRAN (R 4.1.0)
#  tibble         3.1.2      2021-05-16 [1] CRAN (R 4.1.0)
#  tidyr          1.1.3      2021-03-03 [1] CRAN (R 4.1.0)
#  tidyselect     1.1.1      2021-04-30 [1] CRAN (R 4.1.0)
#  usethis      * 2.0.1      2021-02-10 [1] CRAN (R 4.1.0)
#  utf8           1.2.1      2021-03-12 [1] CRAN (R 4.1.0)
#  vctrs          0.3.8      2021-04-29 [1] CRAN (R 4.1.0)
#  withr          2.4.2      2021-04-18 [1] CRAN (R 4.1.0)
#  zip            2.2.0      2021-05-31 [1] CRAN (R 4.1.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
