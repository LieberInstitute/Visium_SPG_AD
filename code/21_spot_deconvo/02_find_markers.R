library("SingleCellExperiment")
library("DeconvoBuddies")
library("tidyverse")
library("sessioninfo")
library("here")
library("HDF5Array")
library("spatialLIBD")
library("cowplot")

sce_in = here('processed-data', '21_spot_deconvo', 'sce_mathys.rds')
spe_in = here('processed-data', '04_build_spe', 'spe_wholegenome.rds')
out_dir = here('processed-data', '21_spot_deconvo')
plot_dir = here('plots', '21_spot_deconvo')

cell_type_var = "broad.cell.type"
find_markers_model = "~individualID"

#   Number of marker genes to use per cell type
n_markers_per_type = 25

################################################################################
#   Functions
################################################################################

write_markers <- function(n_markers, out_path) {
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
        message(num_markers_table)
    }

    stopifnot(all(num_markers_table$num_markers > 0))

    #   Write list of markers
    writeLines(marker_stats_temp$gene, con = out_path)
}

my_plot_expression <- function(
        sce, genes, assay = "logcounts", ct = "cellType", fill_colors = NULL,
        title = NULL, marker_stats
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

    if (!is.null(fill_colors)) {
        expression_violin <- expression_violin + scale_fill_manual(
            values = fill_colors
        )
    }

    # expression_violin
    return(expression_violin)
}

#   Plot mean-ratio distribution by cell type/ layer
boxplot_mean_ratio <- function(n_markers, plot_path, colors_col) {
    p <- marker_stats |>
        filter(rank_ratio <= n_markers) |>
        mutate(ratio, ratio = log(ratio)) |>
        ggplot(aes(cellType.target, ratio, color = cellType.target)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_color_manual(values = metadata(sce)[[colors_col]]) +
        labs(y = "log(Mean Ratio)") +
        theme_bw(base_size = 25) +
        guides(color = "none")

    ggsave(p, filename = plot_path, height = 10, width = 10)
}

################################################################################
#   Find and save markers
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

load(sce_in, verbose = TRUE)

message("Running getMeanRatio2 and findMarkers_1vAll to rank genes as markers...")
marker_stats = get_mean_ratio2(
    sce, cellType_col = cell_type_var, assay_name = "counts"
)
marker_stats_1vall <- findMarkers_1vAll(
    sce, cellType_col = cell_type_var, assay_name = "counts",
    mod = find_markers_model
)
marker_stats <- left_join(
    marker_stats, marker_stats_1vall,
    by = c("gene", "cellType.target")
)

#   Save 'marker_stats' table and the markers themselves (just gene symbols)
saveRDS(marker_stats, file.path(out_dir, 'marker_stats.rds'))
write_markers(marker_stats, file.path(out_dir, 'markers.txt'))

################################################################################
#   Visually check quality of markers
################################################################################

#   Explore how markers look for each cell type
plot_list <- lapply(
    sce[[cell_type_var]],
    function(ct) {
        genes <- marker_stats |>
            filter(
                rank_ratio <= n_markers_per_type,
                cellType.target == ct,
                ratio > 1
            ) |>
            pull(gene)
        my_plotExpression(
            sce, genes,
            ct = cell_type_var,
            fill_colors = metadata(sce)[[colors_col]],
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
