#### load relevant packages ####
library('readxl')
library('dplyr')
library('biomaRt')
library('sessioninfo')
library('here')
library('scran')


# Number of sets: 6 cell types * 4 models = 24 sets * direction (2) = 48 sets
#
# Note: snRNA-seq
#
# Direction available: IndModel.FC or MixedModel.z for the indicator model (Wilcoxon test) or for the mixed effects model (Poisson)
# Statistics available: Yes, but they also took into account effect size with DEGs.Ind.Model and DEGs.Ind.Mix.models
# Ignore direction: 6 cell types * 4 models = 24 sets
# Use direction: 6 cell types * 4 models * 2 directions = 48 sets
# For each of the above, check how many genes we have in each set. They might be too small.
# Note that Excel might introduce some issues in this table too.


input_dir <- "raw-data/GeneSets/2_snRNA-seq/1_Mathys et al_PFC/Mathys et al.xlsx"
mathys_ex <- read_excel(input_dir, sheet = 'Ex')
mathys_in <- read_excel(input_dir, sheet = 'In')
mathys_ast <- read_excel(input_dir, sheet = 'Ast')
mathys_oli <- read_excel(input_dir, sheet = 'Oli')
mathys_opc <- read_excel(input_dir, sheet = 'Opc')
mathys_mic <- read_excel(input_dir, sheet = 'Mic')




