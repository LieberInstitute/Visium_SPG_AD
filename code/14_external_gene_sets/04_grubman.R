library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')

### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensemble_function.R'))



# Number of sets: 6 cell types * concordance (2 options) = 12 sets + direction 6 sets = 18 sets
#
# Note: snRNA-seq
#
# Direction available: Grubman.LogFC (AD vs Control a priori). There's not that many genes.
#
# Statistics available: No. But we have Concordance.
#
#  For each cell type, use all genes (ignore Concordance)
#  For each cell type, filter to Concordance == TRUE
#  For each cell type, use all genes (ignore Concordance), and separate by direction.

input_data <- "raw-data/GeneSets/2_snRNA-seq/2_Grubman et al_Entorhinal/Grubman et al.xlsx"

#comparison with Mathys
table_s3 <-  read_excel(input_data, sheet = "Supplementary Table 3", skip = 10,
                        col_names = TRUE)
#microglia
table_s4a <-read_excel(input_data, sheet = "Supplementary Table 4a", skip = 6,
col_names = TRUE)

#astrocytes
table_s4b <-read_excel(input_data, sheet = "Supplementary Table 4b", skip = 6,
                       col_names = TRUE)

#neurons
table_s4c <-read_excel(input_data, sheet = "Supplementary Table 4c", skip = 6,
                       col_names = TRUE)

#oligos
table_s4d <-read_excel(input_data, sheet = "Supplementary Table 4d", skip = 6,
                       col_names = TRUE)

#opc
table_s4e <-read_excel(input_data, sheet = "Supplementary Table 4e", skip = 6,
                       col_names = TRUE)


