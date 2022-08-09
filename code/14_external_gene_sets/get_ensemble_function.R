##### Convert gene names to ENSEMBL IDs ####
library('biomaRt')
#Check for build
# ensembl <- useEnsembl(biomart = "ensembl")
# biomaRt::searchDatasets( mart = ensembl, pattern = "hsapiens")
# dataset              description    version
# 81 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

hsapiens_genes <- getBM(attributes = c("ensembl_gene_id",
                                       "hgnc_symbol"),
                        mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

hsapiens_genes <- as_tibble(hsapiens_genes)

get_ensemble <- function(table, gene_col){

    table_genes <- hsapiens_genes |> filter(hgnc_symbol %in% table[[gene_col]])
    table <- merge(table, table_genes, by.x = gene_col, by.y = "hgnc_symbol", all.x = TRUE)
    table

}

