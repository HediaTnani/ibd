library(Biobase)
library("sleuth")
library(annotables)
library(biomaRt)
library(dplyr)

wd <- getwd()
# setwd("./..")

# sleuth analysis

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

# so <- sleuth_prep(metadata, target_mapping = ttg, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~sex, 'reduced')
so <- sleuth_fit(so, ~sex + Diagnosis, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_lrt <- filter(sleuth_table_lrt, qval <= 0.05)

pca <- plot_pca(so, color_by = 'Diagnosis')

so$pval_aggregate <- FALSE

so <- sleuth_wt(so, 'DiagnosisUlcerative Colitis', 'full')
so <- sleuth_wt(so, 'sexmale', 'full')
sleuth_table_wt <- sleuth_results(so, 'DiagnosisUlcerative Colitis', test_type='wt', which_model='full', show_all = FALSE)
sleuth_table_wt <- filter(sleuth_table_wt, qval <= 0.05)

volcano <- plot_volcano(so, 'DiagnosisUlcerative Colitis', test_type = 'wt', which_model = 'full', sig_level = 0.0000001, point_alpha = 0.2, sig_color = 'red', highlight = NULL)

so$pval_aggregate <- TRUE

sleuth_live(so)

# figure generation





