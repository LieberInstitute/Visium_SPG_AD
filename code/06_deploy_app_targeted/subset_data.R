library("here")
library("SpatialExperiment")

## Data setup

## Download to my laptop
# scp e:/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/07_spot_qc/spe_targeted_postqc.Rdata processed-data/07_spot_qc/
load(here("processed-data", "07_spot_qc", "spe_targeted_postqc.Rdata"))

## Drop some images
imgData(spe_targeted) <- imgData(spe_targeted)[grepl("lowres", imgData(spe_targeted)$image_id), ]

## Save the final object for the shiny app
save(spe_targeted, file = here("code", "06_deploy_app_targeted", "spe_targeted.Rdata"))
