# Readme.md

- Abeta: contains 8 excel csv files, each of which contains the information about segmented Abeta fragments in a given tissue section that has been stained for the two combinations (comb)f Abeta and genes of interest.
    - Comb1: Abeta, PPP3CA, UCHL1, and SST
    - Comb2: Abeta, IDI1, C3, and NINJ1 
    - Brxxxx: donor identifier number, but here refers to a tissue section from a donor with Brxxxx. Br stands for Brain.
- Comb1_Brxxxx: the output of the HALO FISH-IF module applied to a tissue section from Brxxxx to evaluate the combination of Abeta and 3 genes of interest: 
    - Comb1: Abeta, PPP3CA, UCHL1, and SST
- Comb2_Brxxxx: the output of the HALO FISH-IF module applied to a tissue section from Brxxxx to evaluate the combination of Abeta and 3 genes of interest: 
    - Comb2: Abeta, IDI1, C3, and NINJ1 
    - For Comb2: NINJ1 gene expression was examined separately using the more optimized segmentation setting. So, “NoNINJ1” at the end of each file's title indicates that the NINJ1 dataset was not included, whereas “OnlyNINJ1” indicates that the NINJ1 dataset was the only inclusion.
- The 5 tissue sections are arranged in the following order in both slides.
    - Layer # indicates the annotation made for each AD tissue section to define a region of interest for the downstream HALO FISH-IF analysis: 

|     Br3874    |     Br3854   (Layer2)    |
|---------------|--------------------------|
|               |     Br3873   (Layer3)    |
|               |     Br3880   (Layer4)    |
|               |     Br8549   (Layer5)    |