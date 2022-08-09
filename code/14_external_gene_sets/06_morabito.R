library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')

### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensemble_function.R'))




# Number of sets: 6 cell types * direction (3) * 2 assays (ATAC and RNA) = 36 sets
# Note: snATAC-seq and snRNA-seq
# (low priority since it's ATAC-seq with no assigned genes)
#Table S5: try it with and without direction
# Direction available: avg_logFC
# Statistics available: p_val_adj
# Genes are not assigned. We have the Peak coordinates. We could use GenomicRanges::findOverlaps() to compare against the annotation.
#  Let's try subsetting to Peaks that overlap a single gene first (including the UTRs). If they don't overlap genes or overlap to 2 or more genes, we'll ignore them. But it could be that we lose too many peaks that way.
#     Another option is to find the distance to the closest 5' of a gene. Maybe with a maximum distance filter.
#  Table S1: split by celltype, try it with and without direction
# Direction available: avg_logFC
# Statistics available: p_val_adj (looks filtered already)



###which sheet for table S1?
###which sheet for table S5?
