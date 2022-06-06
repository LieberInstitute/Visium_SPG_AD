
####load relevant libraries####
library('here')
library('readxl')
library('SpatialExperiment')
library('spatialLIBD')
library('ggpubr')
library('dplyr')


#### create output directory ####
dir.create(here::here("plots","13_tge_panel"),
           showWarnings = FALSE,
           recursive =TRUE )
dir_plots <- here::here("plots","13_tge_panel")


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

total_wholegenome <- colSums(counts(spe_wholegenome))
total_tge_wholgenome <- colSums(counts(spe_wholegenome[m_wholegenome, ]))
percent_wholegenome <- total_tge_wholgenome / total * 100


#### match tge genes to spe targeted rowdata####

m_targeted <- match(tge_data$gene_name, rowData(spe_targeted)$gene_name)
m_targeted <-  m_targeted[!is.na(m_targeted)]

total_targeted <- colSums(counts(spe_targeted))
total_tge_targeted<- colSums(counts(spe_targeted[m_targeted, ]))
percent_targeted <- total_tge_targeted/ total * 100

#### add percent columns to colData and create dfs. ####
spe_wholegenome$percent_tge <- percent_wholegenome
spe_targeted$percent_tge <- percent_targeted

#check if row names of both vectors follow same order
#identical(rownames(colData(spe_wholegenome)),
#         rownames(colData(spe_targeted)))
#TRUE

df<- tibble(sample_id =spe_wholegenome$sample_id_short,
                percent_wholegenome =spe_wholegenome$percent_tge,
               percent_targeted =spe_targeted$percent_tge)

df<- df |> group_by(sample_id) |> summarize(across(everything(), mean))


#### create plots and save pdf####
pdf(file.path(dir_plots, paste0("percent_reads_assigned.pdf")), width = 14)

ggpaired(df,
         cond1 ="percent_targeted",
         cond2 = "percent_wholegenome",
         line.color = "gray", line.size = 0.4,
         palette = "npg",
         xlab = "spe type",
         ylab = "mean percent tge reads per sample")


dev.off()



