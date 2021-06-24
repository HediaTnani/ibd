library(Biobase)
library(tximport)
library("sleuth")
library(annotables)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(vegan)

wd <- getwd()
# setwd("./..")

## import data

metadata <- read.csv("../protect/SraRunTable.txt")
metadata <- select(metadata, c('Run', 'Diagnosis', 'sex'))
metadata <- rename(metadata, sample = Run)

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

names(abundances_h5) <- paste0("sample", 1:226)

metadata <- mutate(metadata, path = abundances_h5)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
datasets <- listDatasets(ensembl)

ttg <- read.table("data/tx2gene.txt", sep = "\t", col.names = c("target_id", "ens_gene", "ext_gene"))

## aggregate gene transcripts using tximport

txi.kallisto <- tximport(abundances_h5, type = "kallisto", txOut = FALSE, tx2gene=ttg)

# mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
# ttg <- getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"), mart = mart)
# ttg <- rename(ttg, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# ttg <- select(ttg, c("target_id", "ens_gene", "ext_gene"))

## sleuth analysis

# so <- sleuth_prep(metadata, target_mapping = ttg, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~sex, 'reduced')
so <- sleuth_fit(so, ~sex + Diagnosis, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_table_lrt <- filter(sleuth_table_lrt, qval <= 0.05)

so$pval_aggregate <- FALSE

so <- sleuth_wt(so, 'DiagnosisUlcerative Colitis', 'full')
so <- sleuth_wt(so, 'sexmale', 'full')
sleuth_table_wt <- sleuth_results(so, 'DiagnosisUlcerative Colitis', test_type='wt', which_model='full', show_all = FALSE)
#sleuth_table_wt <- filter(sleuth_table_wt, qval <= 0.05)

so$pval_aggregate <- TRUE

sleuth_live(so)

## figure generation

# volcano plot
volcano_data <-  select(sleuth_table_wt, c("ext_gene","qval", "b"))
volcano_data$log_qval <- -log10(volcano_data$qval)

volcano_data$Legend <- cut(volcano_data$b, breaks = c(-Inf, -4.3, 4.3, Inf), labels = c("-", "0", "+"))
count <- 1

volcano_data$Legend <- as.character(volcano_data$Legend)
for (i in volcano_data$Legend) {
  if (volcano_data$log_qval[count] > 40) {
    volcano_data$Legend[count] <- paste(i, "+", sep="")
  } else {
    volcano_data$Legend[count] <- paste(i, "-", sep="")
  }
  count <- count + 1
}

count <- 1
for (i in volcano_data$Legend) {
  if (i == "++") {
    volcano_data$Legend[count] <- "significant positive"
  } else if (i == "+-") {
    volcano_data$Legend[count] <- "insignificant positive"
  } else if (i == "0+") {
    volcano_data$Legend[count] <- "significant neutral"
  } else if (i == "0-") {
    volcano_data$Legend[count] <- "insignificant neutral"
  } else if (i == "-+") {
    volcano_data$Legend[count] <- "signifcant negative"
  } else if (i == "--") {
    volcano_data$Legend[count] <- "insignificant negative"
  }
  count <- count + 1
}

volcano_data$Legend <- as.factor(volcano_data$Legend)

volcano <- ggplot(volcano_data, aes(x=b, y=-log10(qval), color = Legend))

volcano + 
  geom_point(aes(alpha = 0.5)) +
  geom_text(aes(label = ifelse(Legend == "significant positive", ext_gene, '')), hjust = 0, vjust = 0) +
  geom_vline(xintercept = -4.3, linetype="dashed", size = 0.3) + 
  geom_vline(xintercept = 4.3, linetype="dashed", size = 0.3) + 
  geom_hline(yintercept = 40, linetype="dashed", size = 0.3) + 
  ggtitle("Volcano Plot") + 
  xlab("Log Fold Change") + 
  ylab("Log Odds") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# pca plot



pca <- plot_pca(so, color_by = 'Diagnosis', point_alpha = 0.5)

pca +
  ggtitle("Principal Component Analysis") + 
  xlab("PC1 (62.07%)") + 
  ylab("PC2 (17.01%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# pca variances plot

pca_var <- plot_pc_variance(so, use_filtered = TRUE, units = "est_counts",
                 pca_number = NULL, scale = FALSE, PC_relative = NULL)



