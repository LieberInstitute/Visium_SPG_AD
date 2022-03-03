library("SpatialExperiment")
library("here")
library("spatialLIBD")
library("readxl")
library("RColorBrewer")
library("sessioninfo")

## Load basic SPE data
spe_wholegenome <- readRDS(
    here::here(
        "processed-data", "04_build_spe", "spe_wholegenome.rds"
    )
)
spe_targeted <- readRDS(
    here::here(
        "processed-data", "04_build_spe", "spe_targeted.rds"
    )
)

dir.create(here("plots", "04_build_spe"), showWarnings = FALSE)

## Define order of samples for the grid plots
slide_order <- c("V10A27106", "V10T31036", "V10A27004")
sample_order <- unlist(sapply(slide_order, function(i) {
    sort(unique(spe_wholegenome$sample_id)[grepl(i, unique(spe_wholegenome$sample_id))])
}))
sample_order
#            V10A271061            V10A271062            V10A271063
# "V10A27106_A1_Br3874" "V10A27106_B1_Br3854" "V10A27106_C1_Br3873"
#            V10A271064            V10T310361            V10T310362
# "V10A27106_D1_Br3880" "V10T31036_A1_Br3874" "V10T31036_B1_Br3854"
#            V10T310363            V10T310364            V10A270041
# "V10T31036_C1_Br3873" "V10T31036_D1_Br3880" "V10A27004_A1_Br3874"
#            V10A270042
# "V10A27004_D1_Br3880"


## Edit spatial images to create a black background box
spe_wholegenome <-
    img_update_all(
        spe_wholegenome,
        image_id = "lowres",
        new_image_id = "black",
        brightness = 0,
        saturation = 0,
        hue = 0
    )
imgData(spe_targeted) <- imgData(spe_wholegenome)

## Check the segmentation information
segmentation_variables <-
    c(
        "NAbeta",
        "PAbeta",
        "NDAPI",
        "PDAPI",
        ## These are all 0s right now
        # "NGFAP",
        # "PGFAP",
        # "NLipofuscin",
        # "PLipofuscin",
        # "NMAP2",
        # "PMAP2",
        "NpTau",
        "PpTau"
    )

for (seg_var in segmentation_variables) {
    vis_grid_gene(
        spe_wholegenome,
        geneid = seg_var,
        pdf_file = here(
            "plots",
            "04_build_spe",
            paste0("segmentation_info_", seg_var, ".pdf")
        ),
        spatial = TRUE,
        image_id = "black",
        cont_colors = viridisLite::plasma(101, direction = 1),
        minCount = -1,
        sample_order = sample_order,
        point_size = 2
    )
}

## Make grids for the GraphBased cluster results
length(unique(spe_wholegenome$`10x_graphclust`))
# [1] 9
cols <- RColorBrewer::brewer.pal(length(unique(spe_wholegenome$`10x_graphclust`)), "Set1")
names(cols) <- seq_len(length(cols))
vis_grid_clus(spe_wholegenome,
    clustervar = "10x_graphclust",
    pdf_file = here("plots", "04_build_spe", "wholegenome_graph_based.pdf"),
    sort_clust = TRUE,
    colors = cols,
    spatial = TRUE,
    image_id = "black",
    sample_order = sample_order,
    point_size = 2
)

length(unique(spe_targeted$`10x_graphclust`))
# [1] 6
cols_targeted <- RColorBrewer::brewer.pal(length(unique(spe_targeted$`10x_graphclust`)), "Dark2")
names(cols_targeted) <- seq_len(length(cols_targeted))
vis_grid_clus(spe_targeted,
    clustervar = "10x_graphclust",
    pdf_file = here("plots", "04_build_spe", "targeted_graph_based.pdf"),
    sort_clust = TRUE,
    colors = cols_targeted,
    spatial = TRUE,
    image_id = "black",
    sample_order = sample_order,
    point_size = 2
)

## Read in list of AD genes
ad_genes_raw <- read_xlsx(here("raw-data", "10X_NS targeted gene_AD genes.xlsx"))
ad_genes <- paste0(ad_genes_raw[[2]], "; ", ad_genes_raw[[1]])

addmargins(table(
    "WholeGenome" = ad_genes %in% rowRanges(spe_wholegenome)$gene_search,
    "TargetedSequencing" = ad_genes %in% rowRanges(spe_targeted)$gene_search
))
#            TargetedSequencing
# WholeGenome FALSE TRUE Sum
#       FALSE     3    2   5
#       TRUE      0  134 134
#       Sum       3  136 139
ad_genes[!ad_genes %in% rowRanges(spe_wholegenome)$gene_search]
# [1] "ENSG00000070748; CHAT"    "ENSG00000111537; IFNG"    "ENSG00000113520; IL4"     "ENSG00000254647; INS"
# [5] "ENSG00000187714; SLC18A3"
ad_genes[!ad_genes %in% rowRanges(spe_targeted)$gene_search]
# [1] "ENSG00000070748; CHAT" "ENSG00000111537; IFNG" "ENSG00000113520; IL4"

make_gene_grids <- function(spatial) {
    genes <- ad_genes[ad_genes %in% rowRanges(spatial)$gene_search]
    for (g in genes) {
        print(vis_gene(
            spatial,
            sampleid = "V10A27106_D1_Br3880",
            geneid = g,
            spatial = FALSE,
            assayname = "counts",
            viridis = FALSE
        ))
    }
    return(NULL)
}

pdf(here("plots", "04_build_spe", "wholegenome_V10A27106_D1_Br3880_AD_genes.pdf"), height = 8, width = 9)
make_gene_grids(spe_wholegenome)
dev.off()

pdf(here("plots", "04_build_spe", "targeted_V10A27106_D1_Br3880_AD_genes.pdf"), height = 8, width = 9)
make_gene_grids(spe_targeted)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
