library(Biobase)
library("sleuth")
library(annotables)
library(biomaRt)
library(dplyr)

wd <- getwd()
# setwd("./..")

metadata <- read.csv("../protect/SraRunTable.txt")
metadata <- select(metadata, c('Run', 'Diagnosis', 'sex'))
metadata <- rename(metadata, sample = Run)

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

metadata <- mutate(metadata, path = abundances_h5)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
datasets <- listDatasets(ensembl)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
ttg <- getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"), mart = mart)
ttg <- rename(ttg, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- select(ttg, c("target_id", "ens_gene", "ext_gene"))

so <- sleuth_prep(metadata, target_mapping = ttg, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~sex, 'reduced')

so <- sleuth_fit(so, ~sex + Diagnosis, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- filter(sleuth_table_gene, qval <= 0.05)



# tsv_matrix <- matrix(, nrow = 199240, ncol = 226)

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

# write.csv(tsv_matrix, "tpm_table.csv")

setwd(wd)

