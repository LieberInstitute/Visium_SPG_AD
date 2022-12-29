##### Convert gene names to ENSEMBL IDs ####
library("dplyr")
library("biomaRt")
library("EnsDb.Hsapiens.v79")

get_ensembl <- function(table, gene_col, gene_char) {
    gene_col <- enexpr(gene_col)
    gene_sym_list <- as.data.frame(table |> dplyr::select(!!gene_col)) |> na.omit(!!gene_col)
    gene_sym_list <- c(gene_sym_list[, 1])

    genes_and_IDs <- ensembldb::select(EnsDb.Hsapiens.v79,
        keys = gene_sym_list, keytype = "SYMBOL", columns = c("SYMBOL", "GENEID")
    )
    colnames(genes_and_IDs) <- c("symbol", "gene_ensembl_id")
    table <- merge(table, genes_and_IDs, by.x = gene_char, by.y = "symbol", all.x = TRUE)
    table <- table |> distinct(!!gene_col, .keep_all = TRUE)
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
