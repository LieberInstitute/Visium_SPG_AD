This is a [Samui](https://github.com/chaichontat/samui) app, which displays data from the Visium Spatial Proteogenomics (Visium-SPG)  Alzheimer’s Disease (AD) samples in the [LIBD Visium_SPG_AD project](https://research.libd.org/Visium_SPG_AD/). Here, you can interactively visualize individual spots, along with many associated features.

## Features

You can scroll through these annotations or directly query features by name, through the search box above.

### Pathology Groups

The 7 AD pathology categories. Each spot is given one of the following labels:

* `none`: contains no Aβ and no pTau
* `Ab`: containing only Aβ
* `pTau`: containing only pTau
* `both`: containing both Aβ and pTau
* `n_Ab`: directly adjacent to Aβ spots
* `n_pTau`: directly adjacent to pTau spots
* `n_both`: directly adjacent to spots containing both Aβ and pTau

### Spot Coverage

* `PpTau`: percent of the spot with pTau signal.
* `PAbeta:` percent of the spot with Aβ signal.

### Genes

The logcounts of measured genes at the spot level. Experimentally measured genes were subsetted to only include those with nonzero counts in more than 10% of spots.
