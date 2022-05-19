## import required libraries
library("here")
library("SpatialExperiment")
library("scran")
library("scater")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")

## Load basic SPE data
spe <-
    readRDS(here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "wholegenome",
        paste0("spe_harmony_", "wholegenome", ".rds")
    ))

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "wholegenome",
        "clusters_BayesSpace"
    ),
    prefix = ""
)

## output directories
dir_plots <- here::here("plots", "07_spot_qc", "outliers")
dir.create(dir_plots, showWarnings = FALSE)

path_df <- data.frame(
    spot_id = rownames(colData(spe)),
    diagnosis = colData(spe)$diagnosis,
    sample_id = colData(spe)$sample_id_short,
    NAbeta = as.integer(colData(spe)$NAbeta),
    NpTau = as.integer(colData(spe)$NpTau),
    PAbeta = colData(spe)$PAbeta,
    PpTau = colData(spe)$PpTau
)

n_Abeta <- 1
p_Abeta <- 0.108


path_df <- path_df |> mutate(
    Abeta_outliers = case_when(
        NAbeta > n_Abeta & PAbeta > p_Abeta ~ "n and p",
        NAbeta > n_Abeta & PAbeta <= p_Abeta ~ "n",
        NAbeta <= n_Abeta & PAbeta > p_Abeta ~ "p",
        NAbeta <= n_Abeta & PAbeta <= p_Abeta ~ "none"
    )
)

n_pTau <- 8
p_pTau <- 0.0143

path_df <- path_df |> mutate(
    pTau_outliers = case_when(
        NpTau > n_pTau & PpTau > p_pTau ~ "n and p",
        NpTau > n_pTau & PpTau <= p_pTau ~ "n",
        NpTau <= n_pTau & PpTau > p_pTau ~ "p",
        NpTau <= n_pTau & PpTau <= p_pTau ~ "none"
    )
)

path_df$diagnosis <-
    factor(path_df$diagnosis, levels = c("Control", "AD"))
path_df$Abeta_outliers <- factor(path_df$Abeta_outliers,
                                 levels = c("n", "p", "n and p", "none"))
path_df$pTau_outliers <- factor(path_df$pTau_outliers,
                                levels = c("n", "p", "n and p", "none"))

## Create scatter plots to view outlier information
create_plots <- function(pathology) {
    # pathology can be 'Abeta' or 'pTau'
    colors_hex <- c("#E69F00", "#56B4E9", "#D55E00", "#999999")
    #palette.colors(NULL, "Okabe-Ito")
    #the colorblindr packages suggests using colors from the Okabe-Ito palette
    #for colorblind-friendliness.


    if (pathology == "Abeta") {
        plot <- ggpubr::ggscatter(
            path_df,
            x = "NAbeta",
            y = "PAbeta",
            color = "Abeta_outliers",
            size = 0.5,
            xlab = "Number of Abeta per spot (n)",
            ylab = "Proportion of Abeta per spot (p)"
        ) + guides(colour = guide_legend(override.aes = list(size = 2)))
        plot <- facet(
            plot + theme_bw(base_size = 20) +
                theme(strip.text.x = element_text(size = 20),
                      panel.spacing.x = unit(7, "mm")),

            facet.by = "sample_id",
            short.panel.labs = TRUE
        )
        plot <- plot + scale_color_manual(name = "",
                                          values = colors_hex)
    }

    if (pathology == "pTau") {
        plot <- ggpubr::ggscatter(
            path_df,
            x = "NpTau",
            y = "PpTau",
            color = "pTau_outliers",
            size = 0.5,
            xlab = "Number of pTau per spot (n)",
            ylab = "Proportion of pTau per spot (p)"
        ) + guides(colour = guide_legend(override.aes = list(size = 2)))

        plot <- facet(
            plot + theme_bw(base_size = 20) +
                theme(strip.text.x = element_text(size = 20),
                      panel.spacing.x = unit(7, "mm")),
            facet.by = "sample_id",
            short.panel.labs = TRUE,

        )

        plot <- plot + scale_color_manual(name = "",
                                          values = colors_hex)
    }

    return(plot)
}



# create plots for whole genome
for (path in c("Abeta", "pTau")) {
    pdf(file.path(
        dir_plots,
        paste0(path, "_", "scatter_", "by_", "sample", ".pdf")
    ), width = 14)
    print(create_plots(path))
    dev.off()
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
