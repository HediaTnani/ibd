library(Biobase)
library("sleuth")
library(annotables)
library(biomaRt)
library(dplyr)

wd <- getwd()

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
ttg <- getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"), mart = mart)
ttg <- rename(ttg, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- select(ttg, c("target_id", "ens_gene", "ext_gene"))

this <- read.csv("../protect/tsv/SRR6467548/abundance.tsv", sep = "\t")
transcript_ids = this[,1]
transcript_ids <- sub("[.][0-9]*", "", transcript_ids)
results = getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id'), filters = "ensembl_transcript_id", values = transcript_ids, mart = mart)

metadata <- read.csv("../protect/SraRunTable.txt")

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

tsv_matrix <- matrix(, nrow = 199240, ncol = 226)

patient <- 1
gene <- 1
for (abundance in abundances_tsv) {
  tsv <- read.csv(abundance, sep = "\t")
  while (gene < 199240) {
    tsv_matrix[gene, patient] <- tsv[gene,5]
    gene <- gene + 1
  }
  patient <- patient + 1
  gene <- 1
}
