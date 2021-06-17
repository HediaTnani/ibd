library(Biobase)
library("sleuth")
library(annotables)
library(biomaRt)
library(dplyr)

wd <- getwd()

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

this <- read.csv("../protect/tsv/SRR6467548/abundance.tsv", sep = "\t")
transcript_ids = this[,1]
transcript_ids <- sub("[.][0-9]*", "", transcript_ids)
ttg = getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id'), filters = "ensembl_transcript_id", values = transcript_ids, mart = mart)
gene_ids = ttg$ensembl_gene_id

tte = getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id',"hgnc_symbol","entrezgene_id"), filters = "ensembl_gene_id", values = gene_ids, mart = mart)

different.names <- (!tte$ensembl_transcript_id %in% ttg$ensembl_transcript_id)
not.in.ttg <- tte[different.names]

# metadata <- read.csv("../protect/SraRunTable.txt")
# 
# tsv_folder <- "../protect/tsv"
# tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
# abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
# abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

# tsv_matrix <- matrix(, nrow = 199240, ncol = 226)
# 
# patient <- 1
# gene <- 1
# for (abundance in abundances_tsv) {
#   tsv <- read.csv(abundance, sep = "\t")
#   while (gene < 199240) {
#     tsv_matrix[gene, patient] <- tsv[gene,5]
#     gene <- gene + 1
#   }
#   patient <- patient + 1
#   gene <- 1
# }

rownames(tsv_matrix) <- transcript_ids

# tsv_matrix_rm <- tsv_matrix[grepl("^NA", rownames(df))==F,]

write.csv(tsv_matrix, "../tpm_table_enst.csv")
