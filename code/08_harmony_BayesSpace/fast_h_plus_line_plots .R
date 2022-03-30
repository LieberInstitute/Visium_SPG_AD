library(ggplot2)
library(dplyr)
library(here)


fhplus_data <- read.csv(
    here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        "fasthplus_results.csv"
    ), sep="\t",  #data is tab delimited
)

head(fhplus_data)

fhplus_data<- fhplus_data[fhplus_data$k != "k" ,]

type_list <- c('wholegenome', 'targeted')
spots_set_list <- c('grey_matter', 'all_spots')

##plot output directory
dir_plots <-
    here::here("plots", "08_harmony_BayesSpace")

for (t in type){
    for (s in spots_set){
        pdf(
            file = here::here(
                "plots",
                "08_harmony_BayesSpace",
                paste0(
                    "fast_h_plus_results",
                    t, s,
                    ".pdf"
                )
            ),
            width = 8
        )

        df_subset <- subset(fhplus_data, type == t& spots_set == s)
        plot<- ggplot(df_subset, aes(x=k, y=t_value, group=1)) +
            geom_line()+
            geom_point()

        print(plot)
        dev.off()

    }
}
