library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())

# system("ln -s ../05_deploy_app_wholegenome/sce_pseudo_pathology_wholegenome.rds sce_pseudo_pathology_wholegenome.rds")
rsconnect::deployApp(
    appFiles = c("app.R", "sce_pseudo_pathology_wholegenome.rds", "initial.R"),
    appName = "Visium_SPG_AD_pseudobulk_AD_pathology_wholegenome",
    account = "libd",
    server = "shinyapps.io"
)
