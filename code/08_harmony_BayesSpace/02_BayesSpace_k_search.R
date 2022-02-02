# sgejobs::job_loop(
#     loops = list(spefile = c(
#         "spe_postqc", "spe_targeted_postqc"
#     )),
#     name = "BayesSpace_k_search",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "80G",
#     task_num = 15
# )

## Required libraries
library("getopt")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("RColorBrewer")

## Specify parameters
spec <- matrix(c(
    "spefile", "s", 2, "character", "SPE file name",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Rename from spe_targeted to spe to simplify the code so it can work with
## either
if (opt$spefile == "spe_harmony_targeted") {
    spe <- readRDS(here::here("processed-data", "08_harmony_BayesSpace", "spe_harmony_targeted.rds"))
    suffix <- "targeted"
} else if (opt$spefile == "spe_harmony_wholegenome") {
    spe <- readRDS(here::here("processed-data", "08_harmony_BayesSpace", "spe_harmony_wholegenome.rds"))
    suffix <- "wholegenome"
}

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", suffix)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", suffix)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clustering_results"), showWarnings = FALSE)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(20220127)

spe = spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)

spe$bayesSpace_temp<-spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

sample_ids <- unique(colData(spe)$sample_id)
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots",paste0("vis_clus_bayesSpace_harmony_",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors =  mycolors
  )
  print(my_plot)
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
