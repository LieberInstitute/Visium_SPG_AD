library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')

### load get_ensemble function
source(here('code/14_external_gene_sets/get_ensemble_function.R'))


# Number of sets: 2 + 9 subpopulations * direction = 2 + 18 = 20 sets
# Note: snRNA-seq
# Table S1: try it with and without direction
# Direction available: logFC
# Statistics available: FDR (filtered to 10%)
# Table S2: split by subpopulation, try it with and without direction
# Direction available: logFC
# Statistics available: globalFDR (filtered to 10%)
# Ignore the comparison for now since some of them are very small


#braak comparison's of interest
# 197:371
# 474:552


leung_1 <- read_excel("raw-data/GeneSets/2_snRNA-seq/3_Leung et al/Leung et al.xlsx",
                      sheet = "Supplementary Table 1", col_names = TRUE)



leung_2 <- read_excel("raw-data/GeneSets/2_snRNA-seq/3_Leung et al/Leung et al.xlsx",
                      sheet = "Supplementary Table 2", col_names = TRUE)

