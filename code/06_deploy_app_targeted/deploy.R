library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
load(here("code", "05_deploy_app_wholegenome", ".deploy_info.Rdata"), verbose = TRUE)
rsconnect::setAccountInfo(
    name = deploy_info$name,
    token = deploy_info$token,
    secret = deploy_info$secret
)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "06_deploy_app_targeted"),
    appFiles = c(
        "app.R",
        "spe.Rdata",
        "Visium_IF_AD_modeling_results.Rdata",
        "sce_pseudo_pathology_targeted.rds"
    ),
    appName = "Visium_IF_AD_Kwon2022_TGE",
    account = "libd",
    server = "shinyapps.io"
)
