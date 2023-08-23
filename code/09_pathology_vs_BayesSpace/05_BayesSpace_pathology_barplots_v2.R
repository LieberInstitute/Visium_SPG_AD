## import required libraries
library("here")
library("sessioninfo")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggrepel")

## output directories
dir_plots <-
    here::here(
        "plots",
        "09_pathology_vs_BayesSpace"
    )
dir.create(dir_plots, showWarnings = FALSE)

## Load pathology colors
source(
    here("code", "colors_pathology.R"),
    echo = TRUE,
    max.deparse.length = 500
)

barplots_spe <- function(suffix) {
    ## Load basic SPE data
    spe <- readRDS(here::here(
        "processed-data",
        "07_spot_qc",
        paste0("spe_", suffix, "_postqc.rds")
    ))

    ## Import BayesSpace clusters
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

    ## Import pathology levels
    spe <- cluster_import(
        spe,
        cluster_dir = here::here(
            "processed-data",
            "09_pathology_vs_BayesSpace",
            "pathology_levels"
        ),
        prefix = ""
    )

    ## Convert from character to a factor, so they appear in the order
    ## we want
    spe$path_groups <-
        factor(
            spe$path_groups,
            levels = c(
                "none",
                "Ab",
                "n_Ab",
                "pTau",
                "n_pTau",
                "both",
                "n_both"
            )
        )

    ## Shorten names
    spe$sample_id_shorter <- gsub("Br", "", spe$sample_id_short)

    ## Drop the 3 controls since they are not interesting for this
    cluster_df <-
        as.data.frame(colData(spe[, !grepl("3874$", spe$sample_id)]))
    bayes_cols <-
        cluster_df |> select(matches("BayesSpace_harmony"))

    pdf(file.path(
        dir_plots,
        paste0("BayesSpace_vs_pathology_barplots_", suffix, ".pdf")
    ), width = 24)

    for (i in bayes_cols) {
        plot <- ggplot(cluster_df, aes(x = i, fill = path_groups)) +
            geom_bar(position = "fill", stat = "count") +
            scale_y_continuous(labels = scales::percent) +
            labs(y = "Percentage", x = colnames(i)) +
            geom_text_repel(
                aes(label = stat(count)),
                stat = "count",
                position = "fill",
                size = 3,
                check_overlap = FALSE
            ) +
            scale_fill_manual(values = colors_pathology) +
            facet_grid(. ~ sample_id_shorter) +
            theme_bw(base_size = 20)

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
