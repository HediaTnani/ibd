library(Biobase)
library("sleuth")
library(annotables)
library(tidyverse)

wd <- getwd()
setwd("./..")

this <- read.csv("protect/tsv/SRR6467548/abundance.tsv", sep = "\t")
metadata <- read.csv("protect/SraRunTable.txt")

tsv_matrix <- matrix(, nrow = 199240, ncol = 226)

tsv_folder <- "protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances <- file.path(tsv_subfolders, "abundance.tsv")

patient <- 1
gene <- 1
for (abundance in abundances) {
  tsv <- read.csv(abundance, sep = "\t")
  while (gene < 199240) {
    tsv_matrix[gene, patient] <- tsv[gene,5]
    gene <- gene + 1
  }
  patient <- patient + 1
  gene <- 1
}

#write.csv(tsv_matrix, "tpm_table.csv")

tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

t2g <- tx2gene()

setwd(wd)

