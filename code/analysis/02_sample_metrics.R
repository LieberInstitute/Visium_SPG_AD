## Automatically style the code in this script:
styler::style_file(here::here("code", "analysis", "02_sample_metrics.R"),
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
    here::here("raw-data", "10x_files", "Lieber_Transfer"),
    here::here("raw-data", "10x_files", "Lieber_Transfer_10x_Alignments")
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
    file = here::here("processed-data", "10x_checks", "sample_metrics.Rdata")
)
write.csv(sample_metrics,
    file = here::here("processed-data", "10x_checks", "sample_metrics.csv")
)


dir.create(here("plots", "10x_checks"), showWarnings = FALSE)
pdf(
    here::here("plots", "10x_checks", "spaceranger_metrics_by_number_of_reads.pdf"),
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
    here::here("plots", "10x_checks", "alignment_check_abby_vs_10x.pdf"),
    useDingbats = FALSE,
    width = 10
)
for(i in colnames(sample_metrics)[-c(1, ncol(sample_metrics))]) {
    p <- ggpaired(sample_metrics, id = "Sample.ID", x = "Alignment", y = i, xlab = "Fiducial Alignment", ylab = i, fill = "Alignment", palette = "Dark2", line.color = "gray", ggtheme = theme_pubr(base_size = 30))
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
#  version  R version 4.1.0 Patched (2021-05-18 r80330)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-07-14
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.1.0)
#  cli           3.0.0   2021-06-30 [2] CRAN (R 4.1.0)
#  colorout      1.2-2   2021-05-25 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace    2.0-2   2021-06-24 [2] CRAN (R 4.1.0)
#  crayon        1.4.1   2021-02-08 [2] CRAN (R 4.1.0)
#  DBI           1.1.1   2021-01-15 [2] CRAN (R 4.1.0)
#  digest        0.6.27  2020-10-24 [2] CRAN (R 4.1.0)
#  dplyr         1.0.7   2021-06-18 [2] CRAN (R 4.1.0)
#  ellipsis      0.3.2   2021-04-29 [2] CRAN (R 4.1.0)
#  fansi         0.5.0   2021-05-25 [2] CRAN (R 4.1.0)
#  generics      0.1.0   2020-10-31 [2] CRAN (R 4.1.0)
#  ggplot2     * 3.3.5   2021-06-25 [2] CRAN (R 4.1.0)
#  glue          1.4.2   2020-08-27 [2] CRAN (R 4.1.0)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 4.1.0)
#  here        * 1.0.1   2020-12-13 [1] CRAN (R 4.1.0)
#  htmltools     0.5.1.1 2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets   1.5.3   2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv        1.6.1   2021-05-07 [2] CRAN (R 4.1.0)
#  jsonlite      1.7.2   2020-12-09 [2] CRAN (R 4.1.0)
#  later         1.2.0   2021-04-23 [2] CRAN (R 4.1.0)
#  lattice       0.20-44 2021-05-02 [3] CRAN (R 4.1.0)
#  lifecycle     1.0.0   2021-02-15 [2] CRAN (R 4.1.0)
#  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.1.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.1.0)
#  pillar        1.6.1   2021-05-16 [2] CRAN (R 4.1.0)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.1.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 4.1.0)
#  promises      1.2.0.1 2021-02-11 [2] CRAN (R 4.1.0)
#  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.1.0)
#  R6            2.5.0   2020-10-28 [2] CRAN (R 4.1.0)
#  Rcpp          1.0.7   2021-07-07 [2] CRAN (R 4.1.0)
#  rlang         0.4.11  2021-04-30 [2] CRAN (R 4.1.0)
#  rmote         0.3.4   2021-05-25 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.1.0)
#  scales        1.1.1   2020-05-11 [2] CRAN (R 4.1.0)
#  servr         0.22    2021-04-14 [1] CRAN (R 4.1.0)
#  sessioninfo * 1.1.1   2018-11-05 [2] CRAN (R 4.1.0)
#  tibble        3.1.2   2021-05-16 [2] CRAN (R 4.1.0)
#  tidyselect    1.1.1   2021-04-30 [2] CRAN (R 4.1.0)
#  utf8          1.2.1   2021-03-12 [2] CRAN (R 4.1.0)
#  vctrs         0.3.8   2021-04-29 [2] CRAN (R 4.1.0)
#  withr         2.4.2   2021-04-18 [2] CRAN (R 4.1.0)
#  xfun          0.24    2021-06-15 [2] CRAN (R 4.1.0)
#
# [1] /users/lcollado/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
