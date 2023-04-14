# Readme.md
- FISHIF_Abeta: settings applied to the HALO FISHIF module to analyze Abeta in each tissue section.
- FISHIF_Comb1: settings applied to the HALO FISHIF module to analyze 3 genes of interest from combination 1 in each tissue section.
    - Comb1: PPP3CA, UCHL1, and SST
- FISHIF_Comb2: settings applied to the HALO FISHIF module to analyze 3 genes of interest from combination 2 in each tissue section.
    - Comb2: IDI1, C3, and NINJ1
        - For Comb2: NINJ1 gene expression was examined separately using the more optimized segmentation setting. So, “NoNINJ1” at the end of each file's title indicates that the NINJ1 dataset was not included, whereas “OnlyNINJ1” indicates that the NINJ1 dataset was the only inclusion.
           - The "NoNINJ1" settings still included the less optimized segmentation settings for NINJ1 gene expression, which needs to be ignored for clarity.
           - The "OnlyNINJ1" settings still included those for IDI1 and C3, which can be ignored for clarity.