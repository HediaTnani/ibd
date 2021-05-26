library(Biobase)
library("DESeq2")
library(ggplot2)

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

write.csv(tsv_matrix, "tpm_table.csv")

#dds <- DESeqDataSetFromMatrix(countData=this,colData=metadata,design=~dex,tidy = TRUE)

setwd(wd)

