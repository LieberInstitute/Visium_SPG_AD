Pathology-level `spatialLIBD` documentation
===========================================

This document describes the layer-level portion of the shiny web application made by the  [`spatialLIBD`](https://bioconductor.org/packages/spatialLIBD) Bioconductor package. Here, “layer-level” refers to “pathology-level” data in which we spatially resolved the Visium-SPG transcriptomic data by 7 different pathological categories of Ab and pTau pathology: they are used interchangeably in this context. You can either find the documentation about this package through [Bioconductor](https://bioconductor.org/packages/spatialLIBD) or at the [`spatialLIBD` documentation website](http://lieberinstitute.github.io/spatialLIBD). Below we explain the options common across tabs and each of the tabs at the layer-level data. As explained in the documentation, the layer-level data is the result of pseudo-bulking the spot-level data using the 7 AD pathology categories to compress it, reduce sparsity, and power more analyses.

## Slides and videos

You might find the following slides useful for understanding the features from this part of the web application. Particularly slides 10-12 and 15-22.

<iframe class="speakerdeck-iframe" frameborder="0" src="https://speakerdeck.com/player/dde92cd6dfc04f9589770e074915658f" title="BioTuring_spatialLIBD" allowfullscreen="true" style="border: 0px; background: padding-box padding-box rgba(0, 0, 0, 0.1); margin: 0px; padding: 0px; border-radius: 6px; box-shadow: rgba(0, 0, 0, 0.2) 0px 5px 40px; width: 100%; height: auto; aspect-ratio: 560 / 420;" data-ratio="1.3333333333333333"></iframe>

These slides were part of our 2021-04-27 webinar for BioTuring that you can watch on YouTube:

<iframe width="560" height="315" src="https://www.youtube.com/embed/S8884Kde-1U" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

A recording of an earlier version of this talk is also available on YouTube.

<iframe width="560" height="315" src="https://www.youtube.com/embed/aD2JU-vUv54" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

You might also be interested in this video demonstration of `spatialLIBD` for the [LIBD rstats club](http://research.libd.org/rstatsclub/). Particularly starting at minute 26 with 25 seconds.

<iframe width="560" height="315" src="https://www.youtube.com/embed/LZ2kvCiRVdM?start=1584" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Raw summary

Before the documentation, this tab displays the [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) object that contains the pathology-level data (layer-level in other `spatialLIBD` apps). It's basically useful to know that the data has been loaded and that you can start navigating the app. If you wish to download this data, use the following command.

```{r}
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.12")

## Download sce data
sce_pseudo <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe")
```

Throughout the rest of this document, we'll refer to this object by the name `sce_pseudo`.

This tab also shows the statistical modeling results, described below, that you can access locally and re-shape using the following code.

```{r}
## Reproduce locally with
modeling_results <- fetch_data("Visium_SPG_AD_Visium_wholegenome_modeling_results")
sig_genes <-
        spatialLIBD::sig_genes_extract_all(
            n = nrow(sce_pseudo),
            modeling_results = modeling_results,
            sce_layer = sce_pseudo
        )
```

## Common options

* `Model results`: the statistical modeling results to use. We computed three different types of models:
  1. `enrichment`: one pathological category against all the other 6 pathological categories. Results in t-statistics. `n_Ab`, for example, represents enriched and depleted genes in the `n_Ab` (`next_Ab`) category against all the other 6 pathological categories based on this enrichment modeling method.
  2. `pairwise`: one pathological category against another one. Results in t-statistics with two-sided p-values. SNF8 in `Ab` vs `both`, for example, compares gene expression of SNF8 between the `Ab` (`Abeta`) and `both` (both `Abeta` and `pTau`) categories. 
  3. `anova`: changes among the pathological categories (adjusting for the mean expression) using the data from all 7 pathological categories (`full`) after dropping pathology-associated spots in the white matter layer (`noWM`).

## Reduced dim

In this panel you can visualize the pathology-level data (`sce_pseudo`) across reduced dimensionality representations derived from the gene expression data from the pathology-level pseudo-bulked data. Select which dimensionality reduction method to use with `Reduced Dimension` (PCA, MDS or with the `scater::runPCA` function which provides more stable results and shows the percent of variable explained). Then use `Color by` to choose which variable to color data by, which can be useful to identify groups of pseudo-bulked samples. The options are:

* `age`: age of death of n = 3 donors with AD.
* `APOe`: APOE genotype of n = 3 donors with AD.
* `BCrating`: Braak and CERAD metrics of n = 3 donors with AD.
  - `Br3854`: Braak VI and CERAD Frequent.
  - `Br3873`: Braak V and CERAD Frequent.
  - `Br3880`: Braak VI and CERAD Frequent.
* `Diagnosis`: clinical diagnosis of n = 3 donors with AD.
* `ncells`: number of spots that were combined when pseudo-bulking.
* `path_groups`: 7 pathological categories of AD-related neuropathology.
* `pmi`: post-mortem interval of n = 3 donors with AD.
* `race`: race of n = 3 donors with AD.
* `rin`: RNA integrity number of n = 3 donors with AD.
* `sample_id`: sample identifier of n = 7 samples from n = 3 donors with AD.
* `sex`: sex of n = 3 donors with AD.
* `subject`: donor brain of n = 3 donors with AD.


```{r}
## Reproduce locally with
scater::plotReducedDim(sce_pseudo)
```

## Model boxplots

This tab allows you to make a boxplot of the `logcounts` gene expression from the spatial domain-level data (`sce_pseudo`) for a given `gene`; you can search your gene by typing either the symbol or the Ensembl gene ID. The model result information displayed in the title of the plot is based on which `model results` you selected and whether you are using the short title version or not (controlled by a checkbox). We provide two different color scales you can use: the color blind friendly `viridis` as well as a custom one we used for the `paper`. Through the `Model test` selector, you can choose which particular comparison to display. For example, `Ab` for the enrichment model means that you would display the results of comparing `Ab` against the rest of the 6 pathological categories. `Ab`-`both` for the pairwise model means that you would display the results of the `Ab` category being greater than the `both` category, while `both`-`Ab` is the reverse scenario. Under pairwise, the unused pathological categories are displayed in gray.


Below the plot you can find the subset of the table of results  (`sig_genes` from earlier), sort the table by the different columns, and download it as a CSV if you want. For more details about what each of these columns mean, check the [`spatialLIBD` vignette documentation](http://LieberInstitute.github.io/spatialLIBD/articles/spatialLIBD.html#extract-significant-genes).

```{r}
## Reproduce locally with
spatialLIBD::layer_boxplot()
```

## Gene Set Enrichment

This tab allows you to upload a CSV file that has a particular format as illustrated [in this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/asd_sfari_geneList.csv). This CSV file should contain:

* one column per gene set of interest labeled as column names on the first row,
* no row names, 
* and human Ensembl gene IDs as values in the cells. 

Once you have uploaded a CSV file following this specific format, you can then check if the genes on each of your gene sets are enriched among the statistics from `model results` (`enrichment`, etc) that have a false discovery rate (FDR) adjusted p-value less than `FDR cutoff` (0.1 by default).

Similar to the `Model boxplots` tab, you can interact with the results table or download it.

```{r}
## Reproduce locally with
spatialLIBD::gene_set_enrichment()
spatialLIBD::gene_set_enrichment_plot()
```

## Spatial registration

If you have a single nucleus or single cell RNA-sequencing (snRNA-seq)  (scRNA-seq) dataset, you might group your cells into clusters. Once you do, you could compress the data by pseudo-bulking (like we did to go from `spe` to `sce_pseudo`). You could then compute `enrichment` (`pairwise`, `anova`) statistics for your cell clusters. If you do so, you can then upload a specially formatted CSV file just like the one in [this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer.csv). This file has:

* column names,
* human Ensembl gene IDs as the row names (first column, no name for the column),
* statistics (numeric values) for the cells.

Once you have uploaded a CSV file following this specific format, you can then assess whether the correlation between your statistics and the ones from our spatial domains for the subset of genes (Ensembl ids) present in both. The resulting heatmap and interactive correlation matrix (which again you can interact with and download) can be useful if you are in the process of labeling your sn/scRNA-seq clusters or simply want to compare them against the spatial domain-specific data we have provided. This can also be used for new spatially-resolved transcriptomics datasets.

Finally, you can change the `Maximum correlation` for visualization purposes on the heatmap as it will change the dynamic range for the colors.

```{r}
## Reproduce locally with
spatialLIBD::layer_stat_cor()
spatialLIBD::layer_stat_cor_plot()
```
