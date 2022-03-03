## import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")
library("tidyr")
library("ggpubr")
library("ggplot2")


## Load basic SPE data
spe <- readRDS(
    here::here(
        "processed-data", "07_spot_qc", "spe_wholegenome_postqc.rds"
    )
)



## output directories
dir_plots <- here::here("plots", "07_spot_qc", "pathology_vs_Bayesspace_cluster_barplots")
dir.create(dir_plots, showWarnings = FALSE)


# import cluster info for whole genome
dir_rdata_whole <- here::here("processed-data", "08_harmony_BayesSpace", "wholegenome") # , suffix

cluster_spe <- cluster_import(
    spe,
    cluster_dir = file.path(dir_rdata_whole, "clusters_BayesSpace"),
    prefix = "imported_"
)






colData(cluster_spe)
# imported_BayesSpace_harmony_k15

cluster_df <- as.data.frame(colData(cluster_spe))

cluster_df <- cluster_df |> mutate(across(matches("BayesSpace_harmony"), factor))

cluster_df <- cluster_df |> mutate(
    pTau_outliers =
        ifelse(NpTau > 7 | PpTau > 0.014, "outlier", "normal")
)
cluster_df <- cluster_df |> mutate(
    Abeta_outliers =
        ifelse(NAbeta > 1 | PAbeta > 0.108, "outlier", "normal")
)

addmargins(table(cluster_df$Abeta_outliers))
# normal outlier     Sum
#  36100    2015   38115
2015 / 38115 * 100
# [1] 5.286633
addmargins(table(cluster_df$pTau_outliers))
# normal outlier     Sum
#  27119   10996   38115
10996 / 38115 * 100
# [1] 28.84953

table(cluster_df$NAbeta > 0 | cluster_df$PAbeta > 0)
# FALSE  TRUE
# 35226  2889
2015 / 2889 * 100
# [1] 69.74732

table(cluster_df$NpTau > 0 | cluster_df$PpTau > 0)
# FALSE  TRUE
# 16217 21898
10996 / 21898 * 100
# [1] 50.21463

## pTau observations
# pTau_observations <- cluster_df |> group_by(pTau_outliers, sample_id,i ) |>
#     summarise(n=n()) |> mutate(pTau_outlier_obs = n)

# Stacked bar plots, add labels inside bars

plot <- ggplot(cluster_df, aes(
    fill = pTau_outliers,
    x = imported_BayesSpace_harmony_k14
)) + # , y = PpTau
    geom_bar(position = "fill", stat = "count") +
    geom_text(aes(label = stat(count)), stat = "count", position = "fill")





bayes_cols <- cluster_df |> select(matches("BayesSpace_harmony"))
path_cols <- cluster_df |> select(matches("_outliers"))
##



pdf(file.path(dir_plots, paste0("spe_whole_barplots_pTau_outliers", ".pdf")), width = 14)

for (i in bayes_cols) {
    plot <- ggplot(cluster_df, aes(x = i, fill = pTau_outliers)) +
        geom_bar(position = "fill", stat = "count") +
        scale_y_continuous(labels = scales::percent) +
        labs(y = "Percentage", x = colnames(i)) +
        geom_text(aes(label = stat(count)), stat = "count", position = "fill", size = 2)


    plot <- facet(plot + theme_bw(),
        facet.by = "sample_id",
        short.panel.labs = FALSE
    )

    print(plot)
}
dev.off()

addmargins(table(cluster_df$pTau_outliers[cluster_df$diagnosis != "Control"])) / sum(cluster_df$diagnosis != "Control") * 100
#   normal   outlier       Sum
# 58.32272  41.67728 100.00000
#

pdf(file.path(dir_plots, paste0("spe_whole_pathology_percent_pTau_outliers", ".pdf")), width = 14)
for (i in colnames(bayes_cols)) {
    plot_df <- cluster_df %>%
        group_by(!!sym(i), sample_id) %>%
        summarize(percent_outlier = sum(pTau_outliers == "outlier") / n() * 100, count = sum(pTau_outliers == "outlier"))

    line_df <- plot_df %>%
        group_by(sample_id) %>%
        summarize(mean_line = mean(percent_outlier))

    plot <- ggplot(plot_df, aes(x = !!sym(i), y = percent_outlier)) +
        geom_point() +
        labs(y = "Pathology percent") +
        geom_text(aes(label = count), size = 2, nudge_y = 5)


    plot <- plot + geom_hline(yintercept = 41.67728, color = "red") +
        geom_hline(data = line_df, aes(yintercept = mean_line), color = "grey40")

    plot <- facet(plot + theme_bw(),
        facet.by = "sample_id",
        short.panel.labs = FALSE
    )

    print(plot)
}
dev.off()





pdf(file.path(dir_plots, paste0("spe_whole_barplots_Abeta_outliers", ".pdf")), width = 14)
for (i in bayes_cols) {
    plot <- ggplot(cluster_df, aes(x = i, fill = Abeta_outliers)) +
        geom_bar(position = "fill", stat = "count") +
        scale_y_continuous(labels = scales::percent) +
        labs(y = "Percentage") +
        geom_text(aes(label = stat(count)), stat = "count", position = "fill", size = 2)

    plot <- facet(plot + theme_bw(),
        facet.by = "sample_id",
        short.panel.labs = FALSE
    )


    print(plot)
}
dev.off()

addmargins(table(cluster_df$Abeta_outliers[cluster_df$diagnosis != "Control"])) / sum(cluster_df$diagnosis != "Control") * 100
#    normal    outlier        Sum
# 92.023563   7.976437 100.000000
#


pdf(file.path(dir_plots, paste0("spe_whole_pathology_percent_Abeta_outliers", ".pdf")), width = 14)
for (i in colnames(bayes_cols)) {
    plot_df <- cluster_df %>%
        group_by(!!sym(i), sample_id) %>%
        summarize(percent_outlier = sum(Abeta_outliers == "outlier") / n() * 100, count = sum(Abeta_outliers == "outlier"))

    line_df <- plot_df %>%
        group_by(sample_id) %>%
        summarize(mean_line = mean(percent_outlier))

    plot <- ggplot(plot_df, aes(x = !!sym(i), y = percent_outlier)) +
        geom_point() +
        labs(y = "Pathology percent") +
        geom_text(aes(label = count), size = 2, nudge_y = 5)


    plot <- plot + geom_hline(yintercept = 7.976437, color = "red") +
        geom_hline(data = line_df, aes(yintercept = mean_line), color = "grey40")

    plot <- facet(plot + theme_bw(),
        facet.by = "sample_id",
        short.panel.labs = FALSE
    )


    print(plot)
}
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
