library(Biobase)
library("DESeq2")
library(ggplot2)
library(bigmemory)

wd <- getwd()
setwd("./..")

this <- read.csv("protect/tsv/SRR6467548/abundance.tsv", sep = "\t")
metadata <- read.csv("protect/SraRunTable.txt")

tsv_matrix <- big.matrix(nrow = 199240, ncol = 226, init = 0, backingfile = "tsv_matrix.bin", descriptorfile = "tsv_matrix.desc")

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

#dds <- DESeqDataSetFromMatrix(countData=this,colData=metadata,design=~dex,tidy = TRUE)

setwd(wd)

