library("here")
library("SpatialExperiment")

## Data setup

## Download to my laptop
# scp e:/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/07_spot_qc/spe_postqc.Rdata processed-data/07_spot_qc/
load(here("processed-data", "07_spot_qc", "spe_postqc.Rdata"))

## Drop some images
table(imgData(spe)$image_id)
table(grepl("lowres", imgData(spe)$image_id))
table(imgData(spe)$image_id[grepl("lowres", imgData(spe)$image_id)])
imgData(spe) <- imgData(spe)[grepl("lowres", imgData(spe)$image_id), ]

## Save the final object for the shiny app
save(spe, file = here("code", "05_deploy_app", "spe.Rdata"))
