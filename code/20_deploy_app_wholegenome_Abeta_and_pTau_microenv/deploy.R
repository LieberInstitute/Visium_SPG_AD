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
# options(rsconnect.packrat = TRUE)
## from http://rstudio.github.io/rsconnect/news/index.html#new-features-1-0-0
## The above looks like is only needed to get more informative error messages
## from what I see at https://github.com/rstudio/rsconnect/issues/934#issuecomment-1654104704
## and from what I saw myself. For example with that option set, I got this error:
#### ## End Task Log #########################################################################################################################################################################
#### Error: Unhandled Exception: Child Task 1322356013 failed: Error building image: Error fetching S4Arrays (1.0.4) source. Error downloading package source. Please update your BioConductor packages to the latest version and try again: <BioconductorPackageSource rep
## which I wasn't getting with renv.
rsconnect::deployApp(
    appDir = here("code", "20_deploy_app_wholegenome_Abeta_and_pTau_microenv"),
    appFiles = c(
        "app.R",
        "spe.Rdata",
        "Visium_SPG_AD_modeling_results.Rdata",
        "sce_pseudo_pathology_wholegenome.rds",
        withr::with_dir(here("code", "05_deploy_app_wholegenome"), dir("www", full.names = TRUE))
    ),
    appName = "Visium_SPG_AD_wholegenome_Abeta_and_pTau_microenv",
    account = "libd",
    server = "shinyapps.io"
)
