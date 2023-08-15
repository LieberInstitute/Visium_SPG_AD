library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)
library(basilisk)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(scran)
library(readxl)

mathys_dir = '/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/mathys'
out_dir = here('processed-data', '21_spot_deconvo')
spe_in = here('processed-data', '04_build_spe', 'spe_wholegenome.rds')
num_cores = 4

#   Those for which there are DEGs computed in the Mathys et al paper
cell_types_mathys = c("Ex", "Oli", "In", "Mic", "Opc", "Ast")

dir.create(out_dir, showWarnings = FALSE)

set.seed(20230815)

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
#   Build SingleCellExperiment
#-------------------------------------------------------------------------------

sce = SingleCellExperiment(
    assays = list(counts = dat),
    colData = pheno
)
rownames(sce) = genes

#-------------------------------------------------------------------------------
#   Exclude Mathys et al DEgs as potential markers
#-------------------------------------------------------------------------------

#   The Mathys et al paper includes cell-type-specific diagnosis-associated
#   DEGs. We dont want to use any of these as markers for any cell types,
#   because it could result in confounding gene-expression profiles with
#   diagnosis-related signal. We want markers that work equally well in AD and
#   controls!

#   Form a list of cell-type-specific diagnosis-associated DEGs from Mathys et
#   al manuscript (specifically https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM4_ESM.xlsx )
mathys_deg_list = list()
for (cell_type in cell_types_mathys) {
    mathys_deg_list[[cell_type]] = read_excel(
            file.path(out_dir, 'mathys_supp_tab_2.xlsx'),
            sheet = cell_type,
            skip = 1
        ) |>
        #   Take DEGs that pass criteria listed in their methods
        filter(DEGs.Ind.Model...8, DEGs.Ind.Mix.models...9) |>
        #   Grab just the gene symbols
        pull(...1)
}
mathys_degs = unlist(mathys_deg_list)

#   Take only cell types for which Mathys et al computed DEGs, and exclude
#   those DEGs so they can't be used as markers downstream
sce = sce[!(genes %in% mathys_degs), sce$broad.cell.type %in% cell_types_mathys]

################################################################################
#   Filter and prepare spatial and single-nucleus objects for deconvolution
################################################################################

#-------------------------------------------------------------------------------
#   Compute log normalized counts in SpatialExperiment
#-------------------------------------------------------------------------------

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
spe = readRDS(spe_in)
reducedDims(spe)$spatial = spatialCoords(spe)

message("Running quickCluster()")
Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
spe <-
    computeSumFactors(
        spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(num_cores)
    )
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))

message("Running logNormCounts()")
spe <- logNormCounts(spe)

#-------------------------------------------------------------------------------
#   Filter objects: overlapping, non-mitochondrial genes
#-------------------------------------------------------------------------------

#   Filter out mitochondrial genes (which in single-nucleus data must be
#   technical artifacts, and therefore don't make meaningful markers or training
#   genes for spot deconvolution)
keep = !grepl('^MT-', rownames(sce))
perc_keep <- 100 * (1 - length(which(keep)) / length(keep))
message(
    paste0(
        "Dropped ", perc_keep, "% of total genes when filtering out ",
        "mitochondrial genes"
    )
)
sce = sce[keep, ]

#   Only take the intersection of genes in spatial and single-nucleus objects
shared_symbols = rowData(spe)$gene_name[
    rowData(spe)$gene_name %in% rownames(sce)
]

perc_keep = round(100 * (1 - length(shared_symbols) / nrow(sce)), 1)
message(
    paste0(
        "Dropped ", perc_keep, "% of potential marker genes ",
        "that were not present in the spatial data"
    )
)

sce = sce[shared_symbols, ]
spe = spe[rowData(spe)$gene_name %in% rownames(sce), ]

#   Since genes now line up and rowData(sce) is empty, just take rowData from
#   'spe'
stopifnot(all(rowData(spe)$gene_name == rownames(sce)))
rowData(sce) = rowData(spe)

#   Use EnsemblID for rownames
rownames(sce) = rowData(sce)$gene_id

#-------------------------------------------------------------------------------
#   Save results
#-------------------------------------------------------------------------------

saveRDS(sce, file.path(out_dir, 'sce_mathys.rds'))
saveRDS(spe, file.path(out_dir, 'sce_norm.rds'))
write_anndata(sce, file.path(out_dir, 'adata_mathys.h5ad'))
write_anndata(spe, file.path(out_dir, 'adata_spatial.h5ad'))

session_info()
