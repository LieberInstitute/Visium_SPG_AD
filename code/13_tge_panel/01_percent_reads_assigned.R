
##load relevant libraries
library('here')
library('readxl')


##load tge data
tge_data <- read_xlsx(path =here::here('raw-data',
                                 '10x_Annotated_Human_Neuroscience_Panel.xlsx'),
                      sheet = 'Gene Information')

