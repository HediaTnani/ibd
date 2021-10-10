library(Biobase)
library(tximport)
library(DESeq2)
library("sleuth")
library(annotables)
library(dplyr)
library(ggplot2)
library(vegan)
library(pheatmap)

## import data

metadata <- read.csv("data/SraRunTable.txt")
metadata <- select(metadata, c('Run', 'Diagnosis', 'sex'))
metadata <- rename(metadata, sample = Run)

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

metadata <- mutate(metadata, path = abundances_h5)

ttg <- read.table("data/tx2gene.txt", sep = "\t", col.names = c("target_id", "ens_gene", "ext_gene"))

## tximport for deseq2

names(abundances_h5) <- metadata$sample
txi_kallisto <- tximport(abundances_h5, type = "kallisto", tx2gene = ttg, txOut = TRUE)
txi_sum <- summarizeToGene(txi_kallisto, tx2gene = select(ttg, c("target_id", "ext_gene")))

## deseq2 analysis

sampleTable <- data.frame(condition = factor(metadata$Diagnosis))
rownames(sampleTable) <- colnames(txi_sum$counts) 
dds <- DESeqDataSetFromTximport(txi_sum, sampleTable, ~condition)
dds <- DESeq(dds) 
res <- results(dds, tidy=TRUE)
# gts <- select(ttg, "ens_gene", "ext_gene")
# names(gts) <- c("row", "symbol")
# res_m <- distinct(merge(res, gts, all.x=TRUE))

top_res <- as.data.frame(res$row)
colnames(top_res) <- c("row")
top_res$padj <- -log10(res$padj)
top_res$log <- res$log2FoldChange
top_res <- filter(top_res, padj >= 30 & (log >= 4 | log <= -4))
top_res <- as.data.frame(top_res$row[rev(order(abs(top_res$log)))])
colnames(top_res) <- c("row")
# labels <- as.vector(top_res[seq(1, nrow(top_res), 1), ])
# labels <- as.vector(top_res$row[1:20])
labels <- read.csv("data/protect_results/volcano_labels.csv")
labels <- as.vector(c("DUOXA2","SAA1","DUOX2","MMP10","CXCR1","DEFB4B","SERPINB4","CXCL5","CXCL17","MMP10","TREM1","BMP3","CYP3A4"))


# write.csv(top_res, "data/protect_results/deseq2_differential.csv")

# export deseq2 results

# write.csv(res, "data/protect_results/deseq2_gene.csv")


## figure generation

# deseq2 volcano plot

volcano_data <-  select(res, c("row", "padj", "log2FoldChange"))
volcano_data <- na.omit(volcano_data)
volcano_data$log_qval <- -log10(volcano_data$padj)
#volcano_data <- aggregate(. ~ symbol, volcano_data, sum)

volcano_data$Legend <- cut(volcano_data$log2FoldChange, breaks = c(-Inf, -4, 4, Inf), labels = c("-", "0", "+"))

count <- 1
volcano_data$Legend <- as.character(volcano_data$Legend)
for (i in volcano_data$Legend) {
  if (volcano_data$log_qval[count] > 30) {
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
    volcano_data$Legend[count] <- "significant negative"
  } else if (i == "--") {
    volcano_data$Legend[count] <- "insignificant negative"
  }
  count <- count + 1
}

volcano_data$Legend <- as.factor(volcano_data$Legend)

volcano <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color = Legend))

volcano + 
  geom_point(alpha = 0.5) +
  # geom_text(aes(label = ifelse(Legend == "significant positive" | Legend == "significant negative", row, '')), hjust = 0, vjust = 0) +
  geom_text(aes(label = ifelse(row %in% labels, row, '')), hjust = 0, vjust = 0,show.legend = FALSE,colour = "black") +
  geom_vline(xintercept = -4, linetype="dashed", size = 0.3) + 
  geom_vline(xintercept = 4, linetype="dashed", size = 0.3) + 
  geom_hline(yintercept = 30, linetype="dashed", size = 0.3) + 
  ggtitle("Wald Test") + 
  xlab(expression(paste(log[2], "(Fold Change)"))) + 
  ylab(expression(paste(-log[10], "(Q-value)"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(text = element_text(family = "Helvetica Neue"))+
  theme(legend.position = "none") +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  scale_color_discrete(name = "Significance",labels = c("", "", "","Down in IBD", "", "Up in IBD"))

ggsave(filename = "deseq2_volcano.png", path = "data/figures/", width = 8, height = 6, device='png', dpi=1000)

# deseq2 pca plot

rld <- vst(dds, blind=TRUE)
pca <- plotPCA(rld, intgroup="condition")
pca +
  ggtitle("Principal Component Analysis (DESeq2)") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# deseq2 heatmap

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(
  mat = rld_cor, 
  main = "Correlation Heatmap"
  )


## sleuth analysis

so <- sleuth_prep(metadata, target_mapping = ttg, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5), gene_mode = TRUE)

so <- sleuth_fit(so, ~sex, 'reduced')
so <- sleuth_fit(so, ~sex + Diagnosis, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_table_lrt <- filter(sleuth_table_lrt, qval <= 0.05)

so$pval_aggregate <- FALSE

so <- sleuth_wt(so, 'DiagnosisUlcerative Colitis', 'full')
so <- sleuth_wt(so, 'sexmale', 'full')
sleuth_table_wt <- sleuth_results(so, 'DiagnosisUlcerative Colitis', test_type = 'wt', which_model = 'full', show_all = FALSE)
#sleuth_table_wt <- filter(sleuth_table_wt, qval <= 0.05)

so$pval_aggregate <- TRUE

# sleuth_live(so)

# export sleuth results

# write.csv(sleuth_table_lrt, "data/protect_results/sleuth_lrt.csv")
# write.csv(sleuth_table_wt, "data/protect_results/sleuth_wt.csv")

# sleuth volcano plot

volcano_data <-  select(sleuth_table_wt, c("ext_gene","qval", "b"))
volcano_data$log_qval <- -log10(volcano_data$qval)

volcano_data$Legend <- cut(volcano_data$b, breaks = c(-Inf, -4, 4, Inf), labels = c("-", "0", "+"))
count <- 1

volcano_data$Legend <- as.character(volcano_data$Legend)
for (i in volcano_data$Legend) {
  if (volcano_data$log_qval[count] > 30) {
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
    volcano_data$Legend[count] <- "significant negative"
  } else if (i == "--") {
    volcano_data$Legend[count] <- "insignificant negative"
  }
  count <- count + 1
}

volcano_data$Legend <- as.factor(volcano_data$Legend)

volcano <- ggplot(volcano_data, aes(x=b, y=-log10(qval), color = Legend))

volcano + 
  geom_point(aes(alpha = 0.5)) +
  geom_text(aes(label = ifelse(Legend == "significant positive" | Legend == "significant negative", ext_gene, '')), hjust = 0, vjust = 0) +
  geom_vline(xintercept = -4, linetype="dashed", size = 0.3) + 
  geom_vline(xintercept = 4, linetype="dashed", size = 0.3) + 
  geom_hline(yintercept = 30, linetype="dashed", size = 0.3) + 
  ggtitle("Volcano Plot (Sleuth)") + 
  xlab(expression(paste(log[2], "(Fold Change)"))) + 
  ylab(expression(paste(-log[10], "(Q-value)"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# sleuth pca plot

pca_tpm <- plot_pca(so, units = "tpm", color_by = 'Diagnosis', point_alpha = 0.5)

pca_tpm +
  ggtitle("Principal Component Analysis (TPM)") + 
  xlab("PC1 (62.127%)") + 
  ylab("PC2 (26.475%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_counts <- plot_pca(so, units = "est_counts", color_by = 'Diagnosis', point_alpha = 0.5)

# permanova <- adonis(select(pca_counts$data, c("PC1", "PC2")) ~ pca_counts$data$Diagnosis, method = "manhattan", perm = 999)
# coef <- coefficients(permanova)["pca_counts$data$Diagnosis1",]

pca_counts +
  ggtitle("Principal Component Analysis (Scaled Reads)") + 
  xlab("PC1 (62.27%)") + 
  ylab("PC2 (17.34%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  # stat_ellipse(aes(x = pca_counts$data$PC1, y = pca_counts$data$PC2), type = "t", level = 0.95, linetype = "dashed", size = 0.3)

# pca variances plots

pca_var_tpm <- plot_pc_variance(so, use_filtered = TRUE, units = "tpm",
                 pca_number = NULL, scale = FALSE, PC_relative = NULL)

pca_var_counts <- plot_pc_variance(so, use_filtered = TRUE, units = "est_counts",
                                pca_number = NULL, scale = FALSE, PC_relative = NULL)

