
library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_pseudo_pathology_wholegenome.rds", "initial.R"),
    appName = "Kwon2022_pseudobulk_AD_pathology_wholegenome",
    account = "libd",
    server = "shinyapps.io"
)
