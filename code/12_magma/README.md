README

MAGMA is a self-contained executable and does not need to be installed. 
If the magma file is placed in a directory that is included in the PATH variable, it can be called from anywhere by simply entering 'magma' followed by the desired arguments at the command line. 
Otherwise, the path to the file must be added (eg. './magma' if it is in the current directory.

Documentation and auxiliary files can be found on the MAGMA site at http://ctglab.nl/software/magma
For questions, error reports, suggestions for improvements and feature requests, please email c.a.de.leeuw@vu.nl


Directory Organization for `12_magma`

`01_Jansen_2019` runs MAGMA step 3: gene set analysis of AD-associated genes from 
https://doi.org/10.1038/s41588-018-0311-9

`02_Lancet_2014` runs MAGMA step 3: gene set analysis of FTD-associated genes from 
doi: 10.1016/S1474-4422(14)70065-1.

`03_Nalls_2019` runs step 1 to 3: gene set analysis of PD-associated genes from 
doi: 10.1016/S1474-4422(19)30320-5.

`04_create_genesets` creates Ab and n_Ab related genesets 

`heatmaps` plots MAGMA association test results. 


`fdr_gene_set`, `pvalues_top_50`, `pvalues_top_100` and `pvalues_top_200` are different gene sets we created based on various statistical significance cutoffs. 


