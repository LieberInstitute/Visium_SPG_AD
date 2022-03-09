## Insall paletteer if necessary
if (!requireNamespace("paletteer", quietly = TRUE)) {
    install.packages("paletteer")
}


## Colorblind colors
colors_pathology <- setNames(
    c(
        "grey90",
        paletteer::paletteer_d("dichromat::SteppedSequential_5")[rep(c(6, 18), each = 2) + c(0, 3)],
        paletteer::paletteer_d("beyonce::X7")[4:5]
    )[c(1:3, 6:7, 4:5)],
    c("none", "Ab+", "next_Ab+", "pT+", "next_pT+", "both", "next_both")
)
