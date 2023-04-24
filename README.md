Visium_SPG_AD
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/377886452.svg)](https://zenodo.org/badge/latestdoi/377886452)

## Overview

<img src="http://research.libd.org/Visium_SPG_AD/img/Br3880_D1_pathology.png" align="left" width="300px" />

Welcome to the `Visium_SPG_AD` project! This project involves several
interactive websites, which are publicly accessible for you to browse
and download.

We spatially resolved transcriptomics of the inferior temporal cortex
(ITC) from postmortem human brain tissue to study local gene expression
profiles of brain microenvironments associated with Alzheimer’s disease
(AD) related neuropathology bearing abnormal protein aggregation of
amyloid beta (Abeta; Ab) and phosphorylated Tau (pTau). Using a total of
10 tissue samples from 3 donors with late-stage Alzheimer’s disease (AD)
and 1 age-matched control (Br3874), we generated proteomic-based,
spatially-resolved transcriptome-scale (SRT) maps of the human ITC in
which the immunofluorescence staining of Ab and pTau was integrated into
the complementary SRT data of the identical tissue sections. With Visium
Spatial Proteogenomics
([Visium-SPG](https://www.10xgenomics.com/products/spatial-proteogenomics)),
we utilized a within-subjects design to study gene expression changes
associated with spatially localized brain neuropathology across the 3
donors with AD (Br3854, Br3873, and Br3880). Here, using *BayesSpace* to
identify the gray-matter associated spots in a data-driven manner, we
discovered statistically significant differentially expressed genes in
the core (Ab spots) and peripheral (next_Ab spots) compartments of
Ab-associated microenvironment in the gray matter of ITC. For more
detailed information about study design and experimental results, please
refer to our [manuscript](https://doi.org/10.1101/2023.04.20.537710).
This work was performed by the [Keri
Martinowich](https://www.libd.org/team/keri-martinowich-phd/), [Kristen
Maynard](https://www.libd.org/team/kristen-maynard-phd/), and [Leonardo
Collado-Torres](http://lcolladotor.github.io/) teams at the [Lieber
Institute for Brain Development](libd.org) as well as [Stephanie
Hicks](https://www.stephaniehicks.com/)’s group from [JHBSPH’s
Biostatistics
Department](https://publichealth.jhu.edu/departments/biostatistics).
This project was done in collaboration with [10x
Genomics](https://www.10xgenomics.com/).

This project involves the
[LieberInstitute/Visium_SPG_AD](https://github.com/LieberInstitute/Visium_SPG_AD)
GitHub repository.

If you tweet about this website, the data or the R package please use
the <code>\#Visium_SPG_AD</code> hashtag. You can find previous tweets
that way as shown
<a href="https://twitter.com/search?q=%23Visium_SPG_AD&src=typed_query">here</a>.

Thank you for your interest in our work!

## Study Design

<img src="http://research.libd.org/Visium_SPG_AD/img/study_overview.png" width="1000px" align="left" />

**Spatial transcriptomics combined with immunodetection of Ab and pTau
in the human inferior temporal cortex (ITC)**. (**A**) Schematic of
experimental design using Visium Spatial Proteogenomics (Visium-SPG) to
investigate the impact of Ab and pTau aggregates on the local
microenvironment transcriptome in the post-mortem human brain. Human ITC
blocks were acquired from 3 donors with AD and 1 age-matched
neurotypical control. Tissue blocks were cryosectioned at 10μm to obtain
2-3 replicates per donor and sections were collected onto individual
capture arrays of a Visium spatial gene expression slide, yielding a
total of 3 gene expression experiments. The entire slide (4 tissue
sections) was stained and scanned using multispectral imaging methods to
detect Ab and pTau immunofluorescence (IF) signals as well as
autofluorescence. Following imaging, tissue sections were permeabilized
and subjected to on-slide cDNA synthesis after which libraries were
generated and sequenced. Transcriptomic data was aligned with the
respective IF image data to generate gene expression maps of the local
transcriptome with respect to Ab plaques and pTau elements, including
neurofibrillary tangles. (**B**) High magnification images show Ab
plaques (white triangles) and various neurofibrillary elements such as
tangles (white arrowheads), neuropil threads (red arrowheads), and
neuritic tau plaques (yellow arrowheads). Lipofuscin (cyan) was
identified through spectral unmixing and pixels confounded with this
autofluorescent signal were excluded from analysis, scale bar, 20μm.
(**C**) ITC tissue block from Br3880 (left) and corresponding spotplots
(right) from the Visium data show gene expression of *MOBP* and
*SNAP25*, which demarcates the border between gray matter (GM) and white
matter (WM), scale bar, 1mm. Color scale indicates spot-level gene
expression in logcounts. (**D**) Image processing and quantification of
Ab and pTau per Visium spot. Ab and pTau signals were thresholded in
their single IF channels for segmentation against autofluorescence
background, including lipofuscin. Thresholded Ab and pTau signals were
aligned to the gene expression map of the same tissue section from
Br3880 and quantified as the proportion of number of pixels per Visium
spot, which is visualized in a spotplot, scale bar, 1mm.

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)
- 🔍 [`samui`](https://github.com/chaichontat/samui)

We provide the following interactive websites, organized by software
labeled with emojis:

- 🔭 `spatialLIBD`
  - [Visium_SPG_AD_wholegenome](https://libd.shinyapps.io/Visium_SPG_AD_wholegenome):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    showing the spatially-resolved Visium-SPG data (n = 10 tissue
    sections from 3 AD donors and 1 age-matched control) with
    statistical results comparing the 7 different spot categories
    bearing Ab and pTau pathology within donors with AD. This version
    shows the whole genome gene expression data.
  - [Visium_SPG_AD_TGE](https://libd.shinyapps.io/Visium_SPG_AD_TGE):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    that the same as the previous one, but with the targeted gene
    expression (TGE) panel data.
  - [Visium_SPG_AD_wholegenome_Abeta_microenv](https://libd.shinyapps.io/Visium_SPG_AD_wholegenome_Abeta_microenv):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    similar to `Visium_SPG_AD_wholegenome`, but with the `Ab` and
    `next_Ab` pathologies collapsed into a single one `Ab_env` for
    studying the Abeta micro environment.
- 👀 `iSEE`
  - [Visium_SPG_AD_pseudobulk_AD_pathology_wholegenome](https://libd.shinyapps.io/Visium_SPG_AD_pseudobulk_AD_pathology_wholegenome):
    [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1) website
    showing the pseudo-bulked spatial data of the 7 different
    pathological categories.
- 🔍 `samui`
  - [Visium SPG AD on
    Samui](https://samuibrowser.com/from?url=visium-spg-ad.s3.amazonaws.com/&s=V10A27106_D1_Br3880_AD&s=V10A27004_A1_Br3874_control&s=V10A27004_D1_Br3880_AD&s=V10A27106_A1_Br3874_control&s=V10A27106_B1_Br3854_AD&s=V10A27106_C1_Br3873_AD&s=V10T31036_A1_Br3874_control&s=V10T31036_B1_Br3854_AD&s=V10T31036_C1_Br3873_AD&s=V10T31036_D1_Br3880_AD):
    [`samui`](https://github.com/chaichontat/samui) website that allows
    to zoom in the raw immunostaining image data at the Visium spot
    level.

### Local `spatialLIBD` apps

If you are interested in running the
[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) applications
locally, you can do so thanks to the
[`spatialLIBD::run_app()`](http://research.libd.org/spatialLIBD/reference/run_app.html),
which you can also use with your own data as shown in our [vignette for
publicly available datasets provided by 10x
Genomics](http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html).

``` r
## Run this web application locally with:
spatialLIBD::run_app()
## You will have more control about the length of the session and memory usage.
## See http://research.libd.org/spatialLIBD/reference/run_app.html#examples.
## See also:
## * https://github.com/LieberInstitute/Visium_SPG_AD/tree/master/code/05_deploy_app_wholegenome
## * https://github.com/LieberInstitute/Visium_SPG_AD/tree/master/code/06_deploy_app_targeted
## * https://github.com/LieberInstitute/Visium_SPG_AD/tree/master/code/18_deploy_app_wholegenome_Abeta_microenv
## 
## You could also use spatialLIBD::run_app() to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.
```

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/Visium_SPG_AD/issues](https://github.com/LieberInstitute/Visium_SPG_AD/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Citing our work

Please cite this [manuscript](https://doi.org/10.1101/2023.04.20.537710)
if you use data from this project.

> Influence of Alzheimer’s disease related neuropathology on local
> microenvironment gene expression in the human inferior temporal
> cortex. Sang Ho Kwon, Sowmya Parthiban, Madhavi Tippani, Heena R
> Divecha, Nicholas J Eagles, Jashandeep S Lobana, Stephen R Williams,
> Michelle Mark, Rahul A Bharadwaj, Joel E Kleinman, Thomas M Hyde,
> Stephanie C Page, Stephanie C Hicks, Keri Martinowich, Kristen R
> Maynard, Leonardo Collado-Torres. bioRxiv 2023.04.20.537710; doi:
> <https://doi.org/10.1101/2023.04.20.537710>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article {Kwon2023.04.20.537710,
        author = {Sang Ho Kwon and Sowmya Parthiban and Madhavi Tippani and Heena R Divecha and Nicholas J Eagles and Jashandeep S Lobana and Stephen R Williams and Michelle Mark and Rahul A Bharadwaj and Joel E Kleinman and Thomas M Hyde and Stephanie C Page and Stephanie C Hicks and Keri Martinowich and Kristen R Maynard and Leonardo Collado-Torres},
        title = {Influence of Alzheimer{\textquoteright}s disease related neuropathology on local microenvironment gene expression in the human inferior temporal cortex},
        elocation-id = {2023.04.20.537710},
        year = {2023},
        doi = {10.1101/2023.04.20.537710},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2023/04/20/2023.04.20.537710},
        eprint = {https://www.biorxiv.org/content/early/2023/04/20/2023.04.20.537710.full.pdf},
        journal = {bioRxiv}
    }

### Cite `spatialLIBD`

Below is the citation output from using `citation('spatialLIBD')` in R.
Please run this yourself to check for any updates on how to cite
**spatialLIBD**.

``` r
print(citation("spatialLIBD")[1], bibtex = TRUE)
#> 
#> Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
#> Maynard KR, Collado-Torres L (2022). "spatialLIBD: an R/Bioconductor
#> package to visualize spatially-resolved transcriptomics data." _BMC
#> Genomics_. doi:10.1186/s12864-022-08601-w
#> <https://doi.org/10.1186/s12864-022-08601-w>,
#> <https://doi.org/10.1186/s12864-022-08601-w>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {spatialLIBD: an R/Bioconductor package to visualize spatially-resolved transcriptomics data},
#>     author = {Brenda Pardo and Abby Spangler and Lukas M. Weber and Stephanie C. Hicks and Andrew E. Jaffe and Keri Martinowich and Kristen R. Maynard and Leonardo Collado-Torres},
#>     year = {2022},
#>     journal = {BMC Genomics},
#>     doi = {10.1186/s12864-022-08601-w},
#>     url = {https://doi.org/10.1186/s12864-022-08601-w},
#>   }
```

Please note that the `spatialLIBD` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing the package.

### Cite `VistoSeg`

To cite [`VistoSeg`](http://research.libd.org/VistoSeg/) please use:

> VistoSeg: processing utilities for high-resolution Visium/Visium-IF
> images for spatial transcriptomics data. Madhavi Tippani, Heena R.
> Divecha, Joseph L. Catallini II, Sang Ho Kwon, Lukas M. Weber, Abby
> Spangler, Andrew E. Jaffe, Stephanie C. Hicks, Keri Martinowich,
> Leonardo Collado-Torres, Stephanie C. Page, Kristen R. Maynard bioRxiv
> 2021.08.04.452489; doi: <https://doi.org/10.1101/2021.08.04.452489>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article {Tippani2021.08.04.452489,
        author = {Tippani, Madhavi and Divecha, Heena R. and Catallini, Joseph L. and Kwon, Sang Ho and Weber, Lukas M. and Spangler, Abby and Jaffe, Andrew E. and Hicks, Stephanie C. and Martinowich, Keri and Collado-Torres, Leonardo and Page, Stephanie C. and Maynard, Kristen R.},
        title = {VistoSeg: processing utilities for high-resolution Visium/Visium-IF images for spatial transcriptomics data},
        elocation-id = {2021.08.04.452489},
        year = {2022},
        doi = {10.1101/2021.08.04.452489},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2022/05/13/2021.08.04.452489},
        eprint = {https://www.biorxiv.org/content/early/2022/05/13/2021.08.04.452489.full.pdf},
        journal = {bioRxiv}
    }

## Data Access

We highly value open data sharing and believe that doing so accelerates
science, as was the case between our
[`HumanPilot`](https://doi.org/10.1038/s41593-020-00787-0) and the
external [`BayesSpace`](https://doi.org/10.1038/s41587-021-00935-2)
projects, documented [on this
slide](https://speakerdeck.com/lcolladotor/hca-la-2022?slide=18).

### Processed Data

[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) also allows
you to access the data from this project as ready to use R objects. That
is, a:

- [`SpatialExperiment`](https://doi.org/10.1093/bioinformatics/btac299)
  object for the Visium-SPG samples (n = 10)

You can use the
[`zellkonverter`](https://bioconductor.org/packages/zellkonverter/)
Bioconductor package to convert any of them into Python
[`AnnData`](https://anndata.readthedocs.io/en/latest/) objects. If you
browse our code, you can find examples of such conversions.

If you are unfamiliar with these tools, you might want to check the
[LIBD rstats club](http://research.libd.org/rstatsclub/#.Y4hWlOzMJUM)
(check and search keywords on the
[schedule](https://docs.google.com/spreadsheets/d/1is8dZSd0FZ9Qi1Zvq1uRhm-P1McnJRd_zxdAfCRoMfA/edit?usp=sharing))
videos and resources.

#### Installing spatialLIBD

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialLIBD` from
[Bioconductor](http://bioconductor.org/) with the following code:

``` r
## Install BiocManager in order to install Bioconductor packages properly
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
## Check that you have a valid R/Bioconductor installation
BiocManager::valid()
## Now install spatialLIBD from Bioconductor
## (this version has been tested on macOS, winOS, linux)
BiocManager::install("spatialLIBD")
## If you need the development version from GitHub you can use the following:
# BiocManager::install("LieberInstitute/spatialLIBD")
## Note that this version might include changes that have not been tested
## properly on all operating systems.
```

### R objects

Using `spatialLIBD` you can access the spatialDLPFC transcriptomics data
from the 10x Genomics Visium platform. For example, this is the code you
can use to access the spatially-resolved data. For more details, check
the help file for `fetch_data()`.

``` r
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.12")
## Download the spot-level data
spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")

## This is a SpatialExperiment object
spe
#> class: SpatialExperiment 
#> dim: 27853 38115 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(27853): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
#>   ENSG00000277196
#> rowData names(7): source type ... gene_type gene_search
#> colnames(38115): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1
#>   TTGTTTGTGTAAATTC-1
#> colData names(113): key sample_id ... APOe path_groups_colors
#> reducedDimNames(15): 10x_pca 10x_tsne ... TSNE_perplexity50.HARMONY
#>   TSNE_perplexity80.HARMONY
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
lobstr::obj_size(spe)
#> 2.29 GB

## Remake the logo image
p_pathology <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "path_groups",
    sampleid = "V10A27106_D1_Br3880",
    colors = spe$path_groups_colors[!duplicated(spe$path_groups_colors)],
    spatial = FALSE,
    ... = " Visium SPG AD\nPathology groups -- made with spatialLIBD"
)
p_pathology
```

<a href="https://libd.shinyapps.io/spatialDLPFC_Visium_Sp_pathology"><img src="http://research.libd.org/Visium_SPG_AD/img/Br3880_D1_pathology.png" width="800px" align="center" /></a>

``` r
## Repeat but for BayesSpace at k = 28 (max 28 clusters, can be less)
k_observed <- unique(spe$BayesSpace_harmony_k28)
p_BayesSpace <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_k28",
    sampleid = "V10A27106_D1_Br3880",
    spatial = FALSE,
    colors = setNames(Polychrome::palette36.colors(length(k_observed)), k_observed),
    ... = " Visium SPG AD\nBayesSpace k28 -- made with spatialLIBD"
)
p_BayesSpace
```

<a href="https://libd.shinyapps.io/spatialDLPFC_Visium_Sp_BayesSpace"><img src="http://research.libd.org/Visium_SPG_AD/img/Br3880_D1_BayesSpace.png" width="800px" align="center" /></a>

### Raw data

You can access all the raw data through
[Globus](http://research.libd.org/globus/)
([jhpce#Visium_SPG_AD](http://research.libd.org/globus/jhpce_Visium_SPG_AD/index.html)).
This includes all the input FASTQ files as well as the outputs from
tools such as
[`SpaceRanger`](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger).
The files are organized following the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
project structure.

## Internal

- JHPCE locations:
  `/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD`
- Slack channel:
  [`libd_visium_if_ad_itg`](https://jhu-genomics.slack.com/archives/C01GJ0F731Q).

### Files:

- `code`: R, python, and shell scripts for running various analyses.
- `plots`: plots generated by R analysis scripts in `.pdf` or `.png`
  format
- `processed-data`
  - `Images`: images used for running `SpaceRanger` and other images
  - `spaceranger`: `SpaceRanger` output files
- `raw-data`
  - `10x_files`: files 10x Genomics transferred to us.
  - `FASTQ`: FASTQ files.
  - `Images`: raw images from the scanner in `.tif` format for each
    Visium-SPG slide (around 8GB each). Each slide contains the image
    for four capture areas.

This GitHub repository is organized along the [*R/Bioconductor-powered
Team Data Science* group
guidelines](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk).
It follows the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
structure.

### Other related files

- Reference transcriptome from 10x Genomics:
  `/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/`
