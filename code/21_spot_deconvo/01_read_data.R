library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)
library(basilisk)
library(here)
library(sessioninfo)

mathys_dir = '/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/mathys'
out_dir = here('processed-data', '21_spot_deconvo')
spe_in = here('processed-data', '04_build_spe', 'spe_wholegenome.rds')

dir.create(out_dir, showWarnings = FALSE)

###############################################################################
#  Functions
###############################################################################

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")

                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)

                #  Write AnnData object to disk
                adata$write(filename = filename)

                return()
            },
            env = zellkonverterAnnDataEnv(),
            sce = sce,
            filename = out_path
        )
    )
}

################################################################################
#   Construct SingleCellExperiment from Mathys data
################################################################################

#-------------------------------------------------------------------------------
#   Gather individual components needed to construct the object
#-------------------------------------------------------------------------------

#   Metadata about individuals from which cells were taken
pd = read.csv(
        file.path(mathys_dir, "snRNAseqPFC_BA10_biospecimen_metadata.csv"),
        as.is = TRUE
    ) |>
    as_tibble()

#   Metadata about cells (colData)
pheno = read.delim(
    file.path(mathys_dir, "filtered_column_metadata.txt"), row.names = 1
)
pheno$individualID = pd$individualID[match(pheno$projid, pd$projid)]

#   Counts
dat = readMM(file.path(mathys_dir, "filtered_count_matrix.mtx"))

#   Gene symbols
genes = read.delim(
    file.path(mathys_dir, "filtered_gene_row_names.txt"), header = FALSE,
    as.is = TRUE
) |>
    as_tibble() |>
    pull(V1)

#-------------------------------------------------------------------------------
#   Build SingleCellExperiment and save R and python versions
#-------------------------------------------------------------------------------

sce = SingleCellExperiment(
    assays = list(counts = dat),
    colData = pheno
)
rownames(sce) = genes

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
spe = readRDS(spe_in)
reducedDims(spe)$spatial = spatialCoords(spe)

saveRDS(sce, file.path(out_dir, 'sce_mathys.rds'))
write_anndata(sce, file.path(out_dir, 'adata_mathys.h5ad'))
write_anndata(spe, file.path(out_dir, 'adata_spatial.h5ad'))

session_info()
