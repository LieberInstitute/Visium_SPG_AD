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


## output directories
dir_plots <-
    here::here(
        "plots",
        "09_pathology_vs_BayesSpace",
        "pathology_vs_Bayesspace_cluster_barplots"
    )
dir.create(dir_plots, showWarnings = FALSE)


barplots_spe <- function(suffix) {
    ## Load basic SPE data
    spe <- readRDS(here::here(
        "processed-data",
        "07_spot_qc",
        paste0("spe_", suffix, "_postqc.rds")
    ))

    # import cluster info
    spe <- cluster_import(
        spe,
        cluster_dir = here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            suffix,
            "clusters_BayesSpace"
        ),
        prefix = ""
    )

    cluster_df <- as.data.frame(colData(spe))
    cluster_df$sample_id <- factor(cluster_df$sample_id, levels = unique(cluster_df$sample_id))

    cluster_df <-
        cluster_df |> mutate(across(matches("BayesSpace_harmony"), factor))

    cluster_df <- cluster_df |> mutate(
        pTau_outliers =
            ifelse(NpTau > 8 | PpTau > 0.0143, "outlier", "normal")
    )
    cluster_df <- cluster_df |> mutate(
        Abeta_outliers =
            ifelse(NAbeta > 1 | PAbeta > 0.108, "outlier", "normal")
    )

    print(addmargins(table(cluster_df$Abeta_outliers)))
    # normal outlier     Sum
    #  36100    2015   38115
    2015 / 38115 * 100
    # [1] 5.286633
    print(addmargins(table(cluster_df$pTau_outliers)))
    # normal outlier     Sum
    # 28419    9696   38115
    9696 / 38115 * 100
    # [1] 25.4388

    print(table(cluster_df$NAbeta > 0 | cluster_df$PAbeta > 0))
    # FALSE  TRUE
    # 35226  2889
    2015 / 2889 * 100
    # [1] 69.74732

    print(table(cluster_df$NpTau > 0 | cluster_df$PpTau > 0))
    # FALSE  TRUE
    # 18567 19548
    9696 / 19548 * 100
    # [1] 49.60098

    ## pTau observations
    # pTau_observations <- cluster_df |> group_by(pTau_outliers, sample_id,i ) |>
    #     summarise(n=n()) |> mutate(pTau_outlier_obs = n)

    # Stacked bar plots, add labels inside bars

    # plot <- ggplot(cluster_df, aes(
    #     fill = pTau_outliers,
    #     x = imported_BayesSpace_harmony_k14
    # )) + # , y = PpTau
    #     geom_bar(position = "fill", stat = "count") +
    #     geom_text(aes(label = stat(count)), stat = "count", position = "fill")





    bayes_cols <- cluster_df |> select(matches("BayesSpace_harmony"))
    path_cols <- cluster_df |> select(matches("_outliers"))

    pdf(file.path(
        dir_plots,
        paste0("spe_", suffix, "_barplots_pTau_outliers", ".pdf")
    ), width = 18)

    for (i in bayes_cols) {
        plot <- ggplot(cluster_df, aes(x = i, fill = pTau_outliers)) +
            geom_bar(position = "fill", stat = "count") +
            scale_y_continuous(labels = scales::percent) +
            labs(y = "Percentage", x = colnames(i)) +
            geom_text(
                aes(label = stat(count)),
                stat = "count",
                position = "fill",
                size = 2
            )


        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )

        print(plot)
    }
    dev.off()

    outliers_thres_pTau <- addmargins(table(cluster_df$pTau_outliers[cluster_df$diagnosis != "Control"])) / sum(cluster_df$diagnosis != "Control") * 100
    print(outliers_thres_pTau)
    # normal  outlier      Sum
    # 61.4114  38.5886 100.0000

    pdf(file.path(
        dir_plots,
        paste0("spe_", suffix, "_pathology_percent_pTau_outliers", ".pdf")
    ), width = 18)
    for (i in colnames(bayes_cols)) {
        plot_df <- cluster_df %>%
            group_by(!!sym(i), sample_id) %>%
            summarize(
                percent_outlier = sum(pTau_outliers == "outlier") / n() * 100,
                count = sum(pTau_outliers == "outlier")
            )

        line_df <- plot_df %>%
            group_by(sample_id) %>%
            summarize(mean_line = mean(percent_outlier))

        plot <-
            ggplot(plot_df, aes(x = !!sym(i), y = percent_outlier)) +
            geom_point() +
            labs(y = "Pathology percent") +
            geom_text(aes(label = count),
                size = 2,
                nudge_y = 5
            )


        plot <-
            plot + geom_hline(yintercept = outliers_thres_pTau[2], color = "red") +
            geom_hline(
                data = line_df,
                aes(yintercept = mean_line),
                color = "grey40"
            )

        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )

        print(plot)
    }
    dev.off()


    pdf(file.path(
        dir_plots,
        paste0("spe_", suffix, "_barplots_Abeta_outliers", ".pdf")
    ), width = 18)
    for (i in bayes_cols) {
        plot <- ggplot(cluster_df, aes(x = i, fill = Abeta_outliers)) +
            geom_bar(position = "fill", stat = "count") +
            scale_y_continuous(labels = scales::percent) +
            labs(y = "Percentage") +
            geom_text(
                aes(label = stat(count)),
                stat = "count",
                position = "fill",
                size = 2
            )

        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )


        print(plot)
    }
    dev.off()

    outliers_thres_Abeta <- addmargins(table(cluster_df$Abeta_outliers[cluster_df$diagnosis != "Control"])) / sum(cluster_df$diagnosis != "Control") * 100
    print(outliers_thres_Abeta)
    #    normal    outlier        Sum
    # 92.023563   7.976437 100.000000


    pdf(file.path(
        dir_plots,
        paste0("spe_", suffix, "_pathology_percent_Abeta_outliers", ".pdf")
    ), width = 18)
    for (i in colnames(bayes_cols)) {
        plot_df <- cluster_df %>%
            group_by(!!sym(i), sample_id) %>%
            summarize(
                percent_outlier = sum(Abeta_outliers == "outlier") / n() * 100,
                count = sum(Abeta_outliers == "outlier")
            )

        line_df <- plot_df %>%
            group_by(sample_id) %>%
            summarize(mean_line = mean(percent_outlier))

        plot <-
            ggplot(plot_df, aes(x = !!sym(i), y = percent_outlier)) +
            geom_point() +
            labs(y = "Pathology percent") +
            geom_text(aes(label = count),
                size = 2,
                nudge_y = 5
            )


        plot <-
            plot + geom_hline(yintercept = outliers_thres_Abeta[2], color = "red") +
            geom_hline(
                data = line_df,
                aes(yintercept = mean_line),
                color = "grey40"
            )

        plot <- facet(plot + theme_bw(),
            facet.by = "sample_id",
            short.panel.labs = FALSE
        )


        print(plot)
    }
    dev.off()
}

## Run the function
barplots_spe("wholegenome")
barplots_spe("targeted")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
