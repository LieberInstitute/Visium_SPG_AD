library(ggplot2)
library(dplyr)
library(here)

##locate and read in csv file
fhplus_data <- read.csv(
    here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "fasthplus_results.csv"
    ),
    sep = "\t",
    #data is tab delimited
)

head(fhplus_data)

#remove redundant lines
fhplus_data <- fhplus_data[fhplus_data$k != "k" , ]
#convert k to class integer so it's ordered in the plot
fhplus_data$k <- as.integer(fhplus_data$k)
fhplus_data$fasthplus <- as.numeric(fhplus_data$fasthplus)
fhplus_data$t_value <- as.integer(fhplus_data$t_value)

dim(fhplus_data)
# [1] 99  5
27 * 2 * 2 ## 27 k values * whole/targeted * GM/all spots
# [1] 108

type_list <- c('wholegenome', 'targeted')
spots_set_list <- c('grey_matter', 'all_spots')

with(fhplus_data, tapply(t_value, paste0(type, "_", spots_set), summary))
# $targeted_all_spots
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1408    1905    1905    1885    1905    1905
#
# $targeted_grey_matter
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1680    1713    1713    1712    1714    1714
#
# $wholegenome_all_spots
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1905    1905    1905    1905    1905    1905
#
# $wholegenome_grey_matter
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1664    1667    1668    1667    1668    1668


##plot output directory
dir_plots <-
    here::here("plots", "08_harmony_BayesSpace", "fasthplus")
#dir.create(dir_plots)


##create line plots
for (t in type_list) {
    for (s in spots_set_list) {
        pdf(
            file = here::here(
                "plots",
                "08_harmony_BayesSpace",
                "fasthplus",
                paste0("fasthplus_results_",
                       t, "_", s,
                       ".pdf")
            ),
            width = 8
        )

        df_subset <- subset(fhplus_data, type == t & spots_set == s)
        df_subset <- na.omit(df_subset)  #some fasthplus values were NA
        plot <- ggplot(df_subset, aes(
            x = k,
            y = 1 - fasthplus,
            group = 1
        )) +
            geom_line() +
            geom_point()

        print(plot)
        dev.off()

    }
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
