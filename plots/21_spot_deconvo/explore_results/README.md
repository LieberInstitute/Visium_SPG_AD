Plots are split into directories `white`, `gray`, and `both`, representing visualizations of `AD` patients subsetted to white matter, gray matter, or both (all spots), respectively.

## Plot interpretations

The following plots are generated nearly identically, but swapping the x-axis and fill of the barplots, which represent cell type and pathology group (in either order). An exact interpretation is given below.

- `*/pathology_barplots_inverted.pdf`: One block (pathology group) in one vertical bar (cell type), contains the proportion of counts of that cell type across all spots that belong to that pathology group. In other words, it's the sum of counts for one cell type and pathology group divided by the sum of counts for that cell type and all pathology groups. Note that this computation is done one sample at a time, and the proportions are averaged across all `AD` samples, so that each sample contributes equally even if it contains more spots total than another sample.

- `*/pathology_barplots.pdf`: One block (cell type) in one vertical bar (pathology group) contains the proportion of counts of that cell type across all spots that belong to that pathology group. In other words, it's the sum of counts for that pathology group and one cell type divided by the sum of counts for that pathology group and all cell types. Note that this computation is done one sample at a time, and the proportions are averaged across all `AD` samples, so that each sample contributes equally even if it contains more spots total than another sample.

The boxplots represent similar information to `*/pathology_barplots.pdf`, with the essential difference being how counts are normalized.

- `*/pathology_boxplots.pdf`: A single point in one box of these plots represents, for the particular cell type in question, the average count of that cell type among all spots belonging to the particular pathology group and one sample. Thus each point is one sample. This enables a direct comparison of the *density* of one cell type between different pathology conditions.
