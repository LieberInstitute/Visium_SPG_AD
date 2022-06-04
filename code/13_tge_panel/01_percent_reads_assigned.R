
####load relevant libraries####
library('here')
library('readxl')
library('SpatialExperiment')
library('spatialLIBD')
library('ggplot2')

####load data####
##load tge data
tge_data <- read_xlsx(path =here::here('raw-data',
                                 '10x_Annotated_Human_Neuroscience_Panel.xlsx'),
                      sheet = 'Gene Information')


## load spe data

spe_wholegenome <-
        readRDS(
            here::here(
                "processed-data",
                "08_harmony_BayesSpace",
                "wholegenome",
                paste0("spe_harmony_", 'wholegenome', ".rds")
            )
        )

spe_targeted <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            "targeted",
            paste0("spe_harmony_", 'targeted', ".rds")
        )
    )


#### match tge genes to spe wholegenome rowdata####

m_wholegenome <- match(tge_data$gene_name, rowData(spe_wholegenome)$gene_name)
m_wholegenome <-  m_wholegenome[!is.na(m_wholegenome)]

total <- colSums(counts(spe_wholegenome))
total_tge <- colSums(counts(spe_wholegenome[m_wholegenome, ]))
percent <- total_tge / total * 100


#### match tge genes to spe targeted rowdata####

m_targeted <- match(tge_data$gene_name, rowData(spe_targeted)$gene_name)
m_targeted <-  m_targeted[!is.na(m_targeted)]

total <- colSums(counts(spe_targeted))
total_tge <- colSums(counts(spe_targeted[m_targeted, ]))
percent <- total_tge / total * 100
