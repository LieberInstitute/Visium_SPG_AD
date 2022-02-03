library("spatialLIBD")
library("Polychrome")
library("here")
library("sessioninfo")

plot_SNN_k10 <- function(suffix) {
    spe <- readRDS(
        here::here(
            "processed-data", "08_harmony_BayesSpace", suffix,
            paste0("spe_harmony_", suffix, ".rds")
        )
    )

    cols <- Polychrome::palette36.colors(length(unique(spe$SNN_k10)))
    names(cols) <- unique(spe$SNN_k10)

    vis_grid_clus(
        spe = spe,
        clustervar = "SNN_k10",
        pdf_file = here::here("plots", "08_harmony_BayesSpace", suffix, "graphbased_SNN_k10.pdf"),
        sort_clust = FALSE,
        colors = cols,
        spatial = FALSE,
        point_size = 2
    )

    pdf(here::here("plots", "08_harmony_BayesSpace", suffix, "graphbased_SNN_k10_cut_at.pdf"), height = 24, width = 36)
    for (k in 4:28) {
        k_nice <- sprintf("%02d", k)

        cols <- Polychrome::palette36.colors(k)
        k_var <- paste0("SNN_k10_k", k)
        names(cols) <- unique(colData(spe)[[k_var]])

        p_list <- vis_grid_clus(
            spe = spe,
            clustervar = k_var,
            return_plots = TRUE,
            sort_clust = FALSE,
            colors = cols,
            spatial = FALSE,
            point_size = 2
        )
        print(cowplot::plot_grid(plotlist = p_list))
    }
    dev.off()

    return(NULL)
}

plot_SNN_k10("wholegenome")
plot_SNN_k10("targeted")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
