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

    spe <- cluster_import(spe,
        cluster_dir = here::here("processed-data", "08_harmony_BayesSpace", suffix, "clusters_graphbased")
    )

    cols <- Polychrome::palette36.colors(k)
    names(cols) <- unique(spe$SNN_k10)

    vis_grid_clus(
        spe = spe,
        clustervar = "SNN_k10",
        pdf_file = here::here("plots", "08_harmony_BayesSpace", suffix, paste0("graphbased_SNN_k10", k_nice, ".pdf")),
        sort_clust = FALSE,
        colors = cols,
        spatial = FALSE,
        point_size = 2
    )

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
