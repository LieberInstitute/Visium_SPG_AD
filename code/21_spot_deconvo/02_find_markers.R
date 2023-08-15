library("SingleCellExperiment")
library("DeconvoBuddies")
library("tidyverse")
library("sessioninfo")
library("here")
library("HDF5Array")
library("spatialLIBD")
library("cowplot")
library("colorspace")
library("readxl")

#   Adds the 'spot_plot' function (and more), a wrapper for 'vis_gene' or
#   'vis_clus' with consistent manuscript-appropriate settings
source(here("code", "21_spot_deconvo", "shared_functions.R"))

sce_in = here('processed-data', '21_spot_deconvo', 'sce_mathys.rds')
spe_in = here('processed-data', '04_build_spe', 'spe_wholegenome.rds')
out_dir = here('processed-data', '21_spot_deconvo')
plot_dir = here('plots', '21_spot_deconvo')
sce_in = here('processed-data', '21_spot_deconvo', 'sce_mathys.rds')

cell_type_var = "broad.cell.type"
find_markers_model = "~individualID"
discrete_cell_palette = "Dark 2"
best_looking_sample_id = "V10A27106_D1_Br3880"

#   Number of marker genes to use per cell type
n_markers_per_type = 25

################################################################################
#   Functions
################################################################################

write_markers <- function(marker_stats, n_markers, out_path) {
    #   Take top N marker genes for each cell type
    marker_stats_temp <- marker_stats |>
        filter(
            rank_ratio <= n_markers,
            ratio > 1
        )

    #   Warn if less than the intended number of markers is used for any cell
    #   type
    num_markers_table <- marker_stats_temp |>
        group_by(cellType.target) |>
        summarize(num_markers = n())

    if (any(num_markers_table$num_markers < n_markers)) {
        warning(
            paste(
                "Used less than", n_markers,
                "markers for at least one cell type."
            )
        )
        message("Number of markers per cell type:")
        print(num_markers_table)
    }

    stopifnot(all(num_markers_table$num_markers > 0))

    #   Write list of markers
    writeLines(marker_stats_temp$gene, con = out_path)
}

my_plot_expression <- function(
        sce, genes, assay = "logcounts", ct = "cellType", title = NULL,
        marker_stats
    ) {
    cat_df <- as.data.frame(colData(sce))[, ct, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))

    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    #   Use gene symbols for labels, not Ensembl ID
    symbols <- rowData(sce)$gene_name[match(genes, rownames(sce))]
    names(symbols) <- genes

    #   Add a data frame for adding mean-ratio labels to each gene
    text_df <- marker_stats
    text_df$ratio <- paste0("Mean ratio: ", round(text_df$ratio, 2))
    text_df$Var1 <- factor(text_df$gene, levels = levels(expression_long$Var1))

    expression_violin <- ggplot(
        data = expression_long, aes(x = cat, y = value, fill = cat)
    ) +
        geom_violin(scale = "width") +
        geom_text(
            data = text_df,
            mapping = aes(
                x = length(unique(sce[[ct]])), y = Inf, fill = NULL,
                label = ratio
            ),
            size = 10, hjust = 1, vjust = 1
        ) +
        scale_fill_discrete_qualitative(palette = discrete_cell_palette) +
        facet_wrap(
            ~Var1,
            ncol = 5, scales = "free_y",
            labeller = labeller(Var1 = symbols)
        ) +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw(base_size = 35) +
        theme(
            legend.position = "None", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic")
        ) +
        stat_summary(fun = median, geom = "crossbar", width = 0.3)

    # expression_violin
    return(expression_violin)
}

#   Plot mean-ratio distribution by cell type
boxplot_mean_ratio <- function(marker_stats, n_markers, plot_path) {
    p <- marker_stats |>
        filter(rank_ratio <= n_markers) |>
        mutate(ratio, ratio = log(ratio)) |>
        ggplot(aes(cellType.target, ratio, color = cellType.target)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_color_discrete_qualitative(palette = discrete_cell_palette) +
        labs(y = "log(Mean Ratio)") +
        theme_bw(base_size = 25) +
        guides(color = "none")

    ggsave(p, filename = plot_path, height = 10, width = 10)
}

################################################################################
#   Find and save markers
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

sce = readRDS(sce_in)

#-------------------------------------------------------------------------------
#   Rank genes as potential markers with DeconvoBuddies functions
#-------------------------------------------------------------------------------

message("Running getMeanRatio2 and findMarkers_1vAll to rank genes as markers...")
marker_stats = get_mean_ratio2(
    sce, cellType_col = cell_type_var, assay_name = "logcounts"
)
marker_stats_1vall <- findMarkers_1vAll(
    sce, cellType_col = cell_type_var, assay_name = "logcounts",
    mod = find_markers_model
)
marker_stats <- left_join(
    marker_stats, marker_stats_1vall,
    by = c("gene", "cellType.target")
)
marker_stats$symbol = rowData(sce)$gene_name[
    match(marker_stats$gene, rownames(sce))
]

#   Save 'marker_stats' table and the markers themselves (just Ensembl IDs)
saveRDS(marker_stats, file.path(out_dir, 'marker_stats.rds'))
write_markers(
    marker_stats, n_markers_per_type, file.path(out_dir, 'markers.txt')
)

################################################################################
#   Visually check quality of markers
################################################################################

#   Explore how markers look for each cell type
plot_list <- lapply(
    unique(sce[[cell_type_var]]),
    function(ct) {
        genes <- marker_stats |>
            filter(
                rank_ratio <= n_markers_per_type,
                cellType.target == ct,
                ratio > 1
            ) |>
            pull(gene)
        my_plot_expression(
            sce,
            genes,
            ct = cell_type_var,
            title = paste("Top", length(genes), "for", ct),
            marker_stats = marker_stats |>
                filter(
                    rank_ratio <= n_markers_per_type,
                    cellType.target == ct,
                    ratio > 1
                )
        )
    }
)

#   Write a multi-page PDF with violin plots for each cell group and all
#   markers
pdf(
    file.path(plot_dir, paste0("marker_gene_violin.pdf")),
    width = 35, height = 35
)
print(plot_list)
dev.off()

#   Plot mean ratio against log fold-change for all genes, split by target cell
#   type and colored by whether each gene will be used as a marker
p <- marker_stats |>
    mutate(
        Marker = case_when(
            rank_ratio <= n_markers_per_type ~ paste0(
                "Top-", n_markers_per_type, " Marker"
            ),
            TRUE ~ "Not Marker"
        )
    ) |>
    ggplot(aes(ratio, std.logFC, color = Marker)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~cellType.target, scales = "free_x", ncol = 4) +
    labs(x = "Mean Ratio") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(col = guide_legend(override.aes = list(size = 2)))

pdf(file.path(plot_dir, "mean_ratio_vs_1vall.pdf"), width = 10)
print(p)
dev.off()

#   Plot mean-ratio distibution by cell type
boxplot_mean_ratio(
    marker_stats, n_markers_per_type,
    file.path(plot_dir, "mean_ratio_boxplot.pdf")
)

#-------------------------------------------------------------------------------
#   Make a grid of plots summarizing how sparsely marker genes for each cell
#   type are expressed spatially. Repeat these plots for different numbers of
#   markers per cell type: 15, 25, 50
#-------------------------------------------------------------------------------

spe = readRDS(spe_in)

for (n_markers in c(15, 25, 50)) {
    plot_list <- list()
    i <- 1

    #   Plot proportion of markers having nonzero expression for each cell type
    for (ct in unique(sce[[cell_type_var]])) {
        #   Get markers for this cell type
        markers <- marker_stats |>
            filter(
                cellType.target == ct,
                rank_ratio <= n_markers,
                ratio > 1
            ) |>
            pull(gene)

        for (sample_id in unique(spe$sample_id)) {
            spe_small <- spe[markers, spe$sample_id == sample_id]

            #   For each spot, compute proportion of marker genes with nonzero
            #   expression
            spe_small$prop_nonzero_marker <- colMeans(
                assays(spe_small)$counts > 0
            )

            plot_list[[i]] <- spot_plot(
                spe_small,
                sample_id = sample_id,
                var_name = "prop_nonzero_marker", include_legend = TRUE,
                is_discrete = FALSE, minCount = 0,
                title = paste0(
                    "Prop. markers w/ nonzero exp:\n", ct, " (", sample_id, ")"
                )
            )

            i <- i + 1
        }
    }
    n_sample <- length(unique(spe$sample_id))
    n_rows <- length(unique(marker_stats$cellType.target))

    write_spot_plots(
        plot_list = plot_list, n_col = n_sample, plot_dir = plot_dir,
        file_prefix = paste0("marker_spatial_sparsity_n", n_markers),
        include_individual = FALSE
    )
}

#-------------------------------------------------------------------------------
#   Plot marker sparsity for 15, 20, 50 markers per cell type for all cell types
#   for the best-looking sample
#-------------------------------------------------------------------------------

plot_list <- list()
max_list <- list()
i <- 1

for (n_markers in c(15, 25, 50)) {
    for (cell_type in unique(sce[[cell_type_var]])) {
        #   Get markers for this cell type
        markers <- marker_stats |>
            filter(
                cellType.target == cell_type,
                rank_ratio <= n_markers,
                ratio > 1
            ) |>
            pull(gene)

        spe_small <- spe[markers, spe$sample_id == best_looking_sample_id]

        #   For each spot, compute proportion of marker genes with nonzero
        #   expression
        spe_small$prop_nonzero_marker <- colMeans(
            assays(spe_small)$counts > 0
        )

        max_list[[i]] <- max(spe_small$prop_nonzero_marker)

        plot_list[[i]] <- spot_plot(
            spe_small,
            sample_id = best_looking_sample_id,
            var_name = "prop_nonzero_marker", include_legend = TRUE,
            is_discrete = FALSE, minCount = 0,
            title = sprintf(
                "Prop. %s markers w/ >0 counts (n = %s markers)",
                cell_type, n_markers
            )
        )

        i <- i + 1
    }
}

max_mat <- matrix(
    unlist(max_list),
    ncol = length(unique(sce[[cell_type_var]])), byrow = TRUE
)

#   Now loop back through the plot list (which will be displayed in 2D)
#   and overwrite the scale to go as high as the largest value in the
#   column. This allows for easy comparison between number of markers for
#   the same cell types
for (i_col in 1:length(unique(sce[[cell_type_var]]))) {
    for (i_row in 1:3) {
        index <- (i_row - 1) * length(unique(sce[[cell_type_var]])) + i_col
        upper_limit <- max(max_mat[, i_col])

        plot_list[[index]] <- overwrite_scale(
            plot_list[[index]],
            upper_limit = upper_limit, min_count = 0
        )
    }
}

write_spot_plots(
    plot_list = plot_list, n_col = length(unique(sce[[cell_type_var]])),
    plot_dir = plot_dir, file_prefix = "sparsity_figure",
    include_individual = FALSE
)

session_info()
