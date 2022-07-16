# library(sgejobs)
# sgejobs::job_single(
#         "02_tge_vs_wg_table",
#         create_shell = TRUE,
#         queue = "bluejay",
#         memory = "10G",
#         command = "Rscript 02_tge_vs_wg_table.R",
#         create_logdir = TRUE
#     )


library(sessioninfo)
library(here)
library(dplyr)

### Cloning to Github from RStudio https://www.youtube.com/watch?v=4s5yAxa99cE

#### wg ====

## TODO read file in
# wg_input <-    read.csv(file = )
# here package can be used to input file
# input location Visium_IF_AD/code/05_deploy_app_wholegenome/Visium_IF_AD_wholegenome_model_results_FDR5perc.csv
# input URL on Github #https://github.com/LieberInstitute/Visium_IF_AD/blob/master/code/05_deploy_app_wholegenome/Visium_IF_AD_wholegenome_model_results_FDR5perc.csv


## Alternatively you can also input the file directly from github using
# github link for input
# library RCurl
# x <- getURL
# Y <- read.csv(text = x)
# https://stackoverflow.com/questions/14441729/read-a-csv-from-github-into-r

###


#### tg ====
