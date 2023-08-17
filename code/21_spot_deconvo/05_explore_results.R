library(tidyverse)
library(here)
library(SpatialExperiment)
library(colorspace)
library(spatialLIBD)
library(sessioninfo)

results_in = here('processed-data', '21_spot_deconvo', 'clusters.csv')
plot_dir = here('plots', '21_spot_deconvo', 'explore_results')

discrete_cell_palette = "Dark 2"

################################################################################
#   Functions
################################################################################

#   Given a tibble with columns 'label' (manual layer label), 'deconvo_tool',
#   'cell_type', and 'count', write a set of barplots to PDF under [plot_dir]
#   with name [filename]. 'ylab' give the y-axis label; 'x_var' is the x-axis
#   variable (as a string); 'fill_var' is the fill variable as a string;
#   'fill_scale' is passed to 'scale_fill_manual(values = [fill_scale])';
#   'fill_lab' is the fill label; 'xlab' is the x-axis label
layer_dist_barplot <- function(
        counts_df, out_path, ylab, x_var, fill_var, fill_lab, xlab, fill_palette
    ) {
    p <- ggplot(
        counts_df,
        aes_string(x = x_var, y = "count", fill = fill_var)
    ) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, fill = fill_lab) +
        scale_fill_discrete_qualitative(palette = fill_palette) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    pdf(out_path, width = 10, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Read in data
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

#   Read in spot-deconvolution results and merge with the original
#   SpatialExperiment object
spe = fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
spe = cluster_import(spe, dirname(results_in), prefix = 'c2l_')
spe$c2l_sample = NULL

results = colData(spe) |>
    as_tibble() |>
    pivot_longer(
        starts_with('c2l_'), values_to = 'count', names_to = 'cell_type'
    )

################################################################################
#   Exploratory plots
################################################################################

#-------------------------------------------------------------------------------
#   Barplots that show how cell type proportions vary by pathology group and
#   vice versa
#-------------------------------------------------------------------------------

norm_results = results |>
    group_by(path_groups) |>
    mutate(count = count / sum(count)) |>
    ungroup()

layer_dist_barplot(
    norm_results, out_path = file.path(plot_dir, 'pathology_barplots.pdf'),
    ylab = 'Cell Type Proportion', x_var = 'path_groups',
    fill_var = 'cell_type', xlab = 'Pathology Group', fill_lab = 'Cell Type',
    fill_palette = discrete_cell_palette
)

norm_results = results |>
    group_by(cell_type) |>
    mutate(count = count / sum(count)) |>
    ungroup()

layer_dist_barplot(
    norm_results, out_path = file.path(plot_dir, 'pathology_barplots_inverted.pdf'),
    ylab = 'Pathology Proportion', x_var = 'cell_type',
    fill_var = 'path_groups', xlab = 'Cell Type', fill_lab = 'Pathology Group',
    fill_palette = discrete_cell_palette
)

#-------------------------------------------------------------------------------
#   Boxplots that show how cell type counts vary by pathology group and
#   vice versa
#-------------------------------------------------------------------------------

#   Average counts of each cell type in pathology group
norm_results = results |>
    group_by(path_groups, sample_id, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup()

#   Create plots for each cell type
plot_list = list()
for (cell_type in unique(norm_results$cell_type)) {
    plot_list[[cell_type]] = ggplot(
        norm_results |> filter(cell_type == {{ cell_type }}),
        aes(x = path_groups, y = count)
    ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.05) +
        labs(
            x = "Pathology Group",
            y = "Average Predicted Count",
        ) +
        #   Facet purely for aesthetic purposes: there is only one cell type
        facet_wrap(~cell_type) +
        theme_bw(base_size = 23) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

pdf(file.path(plot_dir, "pathology_boxplots.pdf"), width = 10)
print(plot_list)
dev.off()

session_info()
