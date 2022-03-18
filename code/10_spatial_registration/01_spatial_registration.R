sgejobs::job_single(
    name = "01_spatial_registration",
    create_shell = TRUE,
    queue = "bluejay",
    memory = "10G",
    command = "Rscript 01_spatial_registration.R",
    create_logdir = TRUE

)




library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)


# # ##create directories
# # dir.create(here::here("processed-data","10_spatial_registration"), showWarnings = FALSE)
# dir.create(here::here("processed-data","10_spatial_registration", "pseudo_bulked"), showWarnings = FALSE)
# dir.create(here::here("processed-data", "10_spatial_registration", "dupCor"), showWarnings = FALSE)
# dir.create(here::here("processed-data","10_spatial_registration", "specific_Ts"), showWarnings = FALSE)
#dir.create(here::here("code", "10_spatial_registration"))

k <-  7 # as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

#load post BayesSpace spe object
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            "wholegenome",
            "spe_harmony_wholegenome.rds"
        )
    )

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "wholegenome",
        "clusters_BayesSpace"
    ),
    prefix = ""
)


## Pseudo-bulk for our current BayesSpace cluster results
# sce_pseudo <- aggregateAcrossCells(
#     spe,
#     DataFrame(
#         BayesSpace = colData(spe)[[k]], #REDO
#         sample_id = spe$sample_id
#     )
# )


spe$PseudoSample = paste0(spe$sample_id, ":", colData(spe)[[paste0("BayesSpace_harmony_k", k_nice)]])
cIndexes = splitit(spe$PseudoSample) # gives you the index for each cluster

umiComb <- sapply(cIndexes, function(ii)
    rowSums(assays(spe)$counts[, ii, drop = FALSE])) #

phenoComb = colData(spe)[!duplicated(spe$PseudoSample), ] #creates new colData dataframe with pseudobulked colData
rownames(phenoComb) = phenoComb$PseudoSample #renames rows of new colData frame to be the clusters
phenoComb = phenoComb[colnames(umiComb),]
phenoComb = DataFrame(phenoComb)

sce_pseudo <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(spe)
    ))

#sce_pseudo <- logNormCounts(sce_pseudo) #size factors <0? from aggregateAcrossCells
# range(spe$size_factor)
# [1]  Inf -Inf
# spe$size_factor[spe$size_factor =='NA'] #size factors
# NULL
#

saveRDS(
    sce_pseudo,
    file = here::here(
        "processed-data",
        "10_spatial_registration",
        "pseudo_bulked",
        paste0("sce_pseudobulked_BayesSpace", k_nice, ".RDS")
    )
)

###############################
##### get mean expression  ####
mat <- logcounts(sce_pseudo)

## filter
gIndex = rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter = mat[gIndex,] #subset matrix on just those genes.  want to remove lowly expressed genes.


#####################
## Build a group model

#convert variables to factors
sce_pseudo$spatial.cluster <-
    as.factor(colData(sce_pseudo)[[paste0("BayesSpace_harmony_k", k_nice)]]) #NAs present
sce_pseudo$age <- as.integer(sce_pseudo$age)
sce_pseudo$sex <- as.factor(sce_pseudo$sex)
sce_pseudo$diagnosis <- as.factor(sce_pseudo$diagnosis)
sce_pseudo$subject <- as.factor(sce_pseudo$subject)
#should we be using other variables, like race etc.?


mod <- with(colData(sce_pseudo),
            model.matrix( ~ 0 + spatial.cluster + age + sex + diagnosis))
colnames(mod) <- gsub('cluster', '', colnames(mod)) #not neccesary

## get duplicate correlation
corfit <- duplicateCorrelation(mat_filter, mod,
                               block = sce_pseudo$subject)
saveRDS(
    corfit,
    file = here::here(
        "processed-data",
        "10_spatial_registration",
        "dupCor",
        paste0("pseudobulked_dupCor_k", k_nice, ".RDS")
    )
)


## Next for each layer test that layer vs the rest
cell_idx <-
    splitit(sce_pseudo$spatial.cluster)  #cell_idx is actually cluster_idx

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudo))
    res[x] <-
        1   #indicator of whether pseudobulked column belongs to
    m <-
        with(colData(sce_pseudo),
             #find genes that are diff expressed across diff BayesSpace clusters adjusting for
             model.matrix( ~ res + diagnosis + age + sex))      #age, diagnosis, sex
    eBayes(
        lmFit(
            mat_filter,
            design = m,
            block = sce_pseudo$sample_id,
            ##block helps take into account that there are some bulks that come from the same sample and
            ##they may be highly correlated to each other
            correlation = corfit$consensus.correlation
        )
    )
})

saveRDS(
    eb0_list_cell,
    file = here::here(
        "processed-data",
        "10_spatial_registration",
        "specific_Ts",
        paste0("pseudobulked_specific_Ts_k", k_nice, ".RDS")
    )
)


##########
## Extract the p-values

pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat_filter)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat_filter)
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
    'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
                           t0_contrasts_cell > 0),
    'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
                                t0_contrasts_cell > 0),
    'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
                                t0_contrasts_cell > 0)
)

#for k = 04
# FDRsig Pval10.6sig Pval10.8sig
# 1  16576        2113         401
# 2    266          12           1
# 3     38           0           0
# 4      0           0           0

###################
load(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata"
)
load(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb0_list.Rdata"
)

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

############
# line up ##
#find genes that match b/w both datasets
mm = match(rownames(pvals0_contrasts), rownames(pvals0_contrasts_cell))

#subset data from Maynard 2021
pvals0_contrasts = pvals0_contrasts[!is.na(mm),]
t0_contrasts = t0_contrasts[!is.na(mm),]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm),]

#subset from current data
pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)],]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)],]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)],]

#correlation b/w t-stats
cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left -Top 100 genes from Maynard 2021
layer_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(t0_contrasts),
as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

#corr matrix for t-stats with Top 100 x 7 layers
cor_t_layer = cor(t0_contrasts_cell[layer_ind,],
                  t0_contrasts[layer_ind,])
signif(cor_t_layer, 3)

### heatmap
theSeq = seq(-.85, .85, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

dd = dist(1 - cor_t_layer)
hc = hclust(dd)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]
colnames(cor_t_layer_toPlot) = gsub("ayer", "", colnames(cor_t_layer_toPlot))
#rownames(cor_t_layer_toPlot)[rownames(cor_t_layer_toPlot) == "Oligodendrocytes"] = "OLIGO" # doesn't matter

##plot output directory
dir_plots <-
    here::here("plots", "10_spatial_registration")
dir.create(dir_plots, showWarnings = FALSE)

#http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html newer function for plotting

pdf(
    file = here::here(
        "plots",
        "10_spatial_registration",
        paste0(
            "pseudobulked_bayesSpace_vs_mannual_annotations_k",
            k_nice,
            ".pdf"
        )
    ),
    width = 8
)
print(
    levelplot(
        cor_t_layer_toPlot,
        aspect = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "",
        scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
    )
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
