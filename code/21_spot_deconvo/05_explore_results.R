library(tidyverse)
library(here)
library(SpatialExperiment)
library(colorspace)
library(spatialLIBD)
library(sessioninfo)
library(getopt)

spec = matrix(
    c("subset", "s", 1, "character", "Tissue type to subset to"),
    byrow = TRUE, ncol = 5
)
opt = getopt(spec)

results_in = here('processed-data', '21_spot_deconvo', 'clusters.csv')
plot_dir = here('plots', '21_spot_deconvo', 'explore_results', opt$subset)

################################################################################
#   Functions
################################################################################

#   Given a tibble with 2 discrete columns and a continuous 'count' column
#   representing spot deconvolution results, write a set of barplots to
#   [out_path]. 'ylab' give the y-axis label; 'x_var' is the x-axis
#   variable (as a string); 'fill_var' is the fill variable as a string;
#   'fill_palette' is passed to
#   'scale_fill_discrete_qualitative(palette = fill_palette)'; 'fill_lab' is the
#   fill label; 'xlab' is the x-axis label. Returns NULL
layer_dist_barplot <- function(
        counts_df, out_path, ylab, x_var, fill_var, fill_lab, xlab, fill_values
    ) {
    p <- ggplot(
        counts_df,
        aes_string(x = x_var, y = "count", fill = fill_var)
    ) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, fill = fill_lab) +
        scale_fill_manual(values = fill_values) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    pdf(out_path, width = 10, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Read in data
################################################################################

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(08212023)

#   Read in spot-deconvolution results and merge with the original
#   SpatialExperiment object. Subset to only AD samples
spe = fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")
spe = cluster_import(spe, dirname(results_in), prefix = 'c2l_')
spe$c2l_sample = NULL
spe = spe[, spe$diagnosis == "AD"]

#   Subset to white matter, gray matter, or neither
if (opt$subset == "white") {
    spe = spe[, spe$BayesSpace_harmony_k02 == 2]
} else if (opt$subset == "gray") {
    spe = spe[, spe$BayesSpace_harmony_k02 == 1]
}

results = colData(spe) |>
    as_tibble() |>
    #   Each row becomes a unique cell type
    pivot_longer(
        starts_with('c2l_'), values_to = 'count', names_to = 'cell_type'
    ) |>
    #   Clean up cell-type names for plotting
    mutate(
        cell_type = case_when(
            cell_type == 'c2l_ast' ~ 'Astro',
            cell_type == 'c2l_ex' ~ 'Excit',
            cell_type == 'c2l_in' ~ 'Inhib',
            cell_type == 'c2l_mic' ~ 'Micro',
            cell_type == 'c2l_oli' ~ 'Oligo',
            cell_type == 'c2l_opc' ~ 'OPC'
        )
    )

#   Verify that input (VistoSeg) and output (cell2location) cell counts per spot
#   are similar
message(sprintf("For '%s' subset and AD samples:", opt$subset))
message(
    sprintf("Mean VistoSeg cell counts per spot: %s", round(mean(spe$NDAPI), 2))
)
c2l_av = results |>
    group_by(key) |>
    summarize(count = sum(count)) |>
    summarize(av = mean(count)) |>
    pull(av) |>
    round(2)
message(sprintf("Mean Cell2location cell counts per spot: %s", c2l_av))

################################################################################
#   Exploratory plots
################################################################################

path_colors = unique(results$path_groups_colors)
names(path_colors) = sapply(
    path_colors,
    function(x) results[match(x, results$path_groups_colors), ][['path_groups']]
)

cell_colors = qualitative_hcl(
    length(unique(results$cell_type)), palette = "Dark 3"
)

#-------------------------------------------------------------------------------
#   Barplots that show how cell type proportions vary by pathology group and
#   vice versa
#-------------------------------------------------------------------------------

norm_results = results |>
    #   For each pathology group and sample_id, normalize by the total counts of
    #   all cell types
    group_by(path_groups, sample_id) |>
    mutate(count = count / sum(count)) |>
    #   Now for each path group, sample_id and cell type, add up counts for all
    #   relevant spots
    group_by(path_groups, cell_type, sample_id) |>
    summarize(count = sum(count)) |>
    #   Now average across samples
    group_by(path_groups, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup()

layer_dist_barplot(
    norm_results, out_path = file.path(plot_dir, 'pathology_barplots.pdf'),
    ylab = 'Cell Type Proportion', x_var = 'path_groups',
    fill_var = 'cell_type', xlab = 'Pathology Group', fill_lab = 'Cell Type',
    fill_values = cell_colors
)

norm_results = results |>
    #   For each cell type and sample_id, normalize by the total counts of
    #   cells in all pathology groups
    group_by(cell_type, sample_id) |>
    mutate(count = count / sum(count)) |>
    #   Now for each path group, sample_id and cell type, add up counts for all
    #   relevant spots
    group_by(path_groups, cell_type, sample_id) |>
    summarize(count = sum(count)) |>
    #   Now average across samples
    group_by(path_groups, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup()

layer_dist_barplot(
    norm_results, out_path = file.path(plot_dir, 'pathology_barplots_inverted.pdf'),
    ylab = 'Pathology Proportion', x_var = 'cell_type',
    fill_var = 'path_groups', xlab = 'Cell Type', fill_lab = 'Pathology Group',
    fill_values = path_colors
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
        aes(x = path_groups, y = count, fill = path_groups)
    ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.05) +
        labs(
            x = "Pathology Group",
            y = "Average Predicted Count",
        ) +
        scale_fill_manual(values = path_colors, guide = "none") +
        #   Facet purely for aesthetic purposes: there is only one cell type
        facet_wrap(~cell_type) +
        theme_bw(base_size = 23) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

pdf(file.path(plot_dir, "pathology_boxplots.pdf"), width = 10)
print(plot_list)
dev.off()

#   For gray matter, check if Excit counts significantly differ in pTau against
#   all other groups combined
if (opt$subset == "gray") {
    stat_results = norm_results |>
        filter(cell_type == 'Excit') |>
        mutate(is_pTau = path_groups == 'pTau')

    t_results = t.test(count ~ is_pTau, stat_results, var.equal = TRUE)
    message(
        paste(
            'p-value for Excit counts in pTau vs. all:',
            signif(t_results$p.value, 3)
        )
    )
}

#-------------------------------------------------------------------------------
#   Boxplots to assess how glial / neuronal ratio changes by pathology group
#-------------------------------------------------------------------------------

norm_results = results |>
    group_by(sample_id, path_groups, cell_type) |>
    summarize(count = mean(count)) |>
    pivot_wider(names_from = cell_type, values_from = count) |>
    group_by(sample_id, path_groups) |>
    summarize(
        glial_neuron = (Astro + Micro) / (Excit + Inhib),
        glial_excit = (Astro + Micro) / Excit
    ) |>
    ungroup()

#   Glial vs. all neurons
p = ggplot(
        norm_results, aes(x = path_groups, y = glial_neuron, fill = path_groups)
    ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.05) +
    labs(x = "Pathology Group", y = "(Astro + Micro) / (Excit + Inhib)") +
    scale_fill_manual(values = path_colors, guide = "none") +
    theme_bw(base_size = 23) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

pdf(file.path(plot_dir, "pathology_glial_v_neuron.pdf"), width = 10)
print(p)
dev.off()

#   Glial vs. excitatory neurons
p = ggplot(
        norm_results, aes(x = path_groups, y = glial_excit, fill = path_groups)
    ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.05) +
    labs(x = "Pathology Group", y = "(Astro + Micro) / Excit") +
    scale_fill_manual(values = path_colors, guide = "none") +
    theme_bw(base_size = 23) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

pdf(file.path(plot_dir, "pathology_glial_v_excit.pdf"), width = 10)
print(p)
dev.off()

#   Test if glial/neuronal ratio significantly differs between n_Ab and all
#   other groups combined
if (opt$subset == "gray") {
    norm_results$is_n_Ab = norm_results$path_groups == 'n_Ab'
    t_results = t.test(glial_neuron ~ is_n_Ab, norm_results, var.equal = TRUE)
    message(
        paste(
            'p-value for (Astro + Micro) / (Excit + Inhib) ratios in n_Ab vs. all:',
            signif(t_results$p.value, 3)
        )
    )
}

#   ANOVA for glial/neuronal ratios among pathology groups
if ((opt$subset == "gray") || (opt$subset == "white")) {
    anov = oneway.test(
        glial_neuron ~ path_groups, norm_results, var.equal = TRUE
    )
    message(
        paste(
            'ANOVA for (Astro + Micro) / (Excit + Inhib) ratios:',
            signif(anov$p.value, 3)
        )
    )
}

session_info()
