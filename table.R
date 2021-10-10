library(Biobase)
library(annotables)
library(dplyr)
library(vegan)
library(data.table)

## produce transcript-level tpm table

ttg <- read.table("data/tx2gene.txt", sep = "\t", col.names = c("target_id", "ens_gene", "ext_gene"))

metadata <- read.csv("data/SraRunTable.txt")

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

row <- read.csv(abundances_tsv[[1]], sep = "\t")$target_id
col <- metadata$Run
df_row <- data.frame(row)
ttg_name <- select(ttg, c("target_id", "ext_gene"))
ttg_id <- select(ttg, c("target_id", "ens_gene"))

names(ttg_name) <- c("row", "symbol")
# ttg_name$symbol[[187688]] <- "5-nucleotidase"
# ttg_name$symbol[[187461]] <- "Wilms"
reduced_gene <- merge(df_row, ttg_name, all.x=TRUE)  

names(ttg_id) <- c("row", "id")
reduced_id <- merge(df_row, ttg_id, all.x=TRUE)  

tsv_df <- data.frame(matrix(ncol=226, nrow=199240, dimnames=list(row, col)))

## transcript-level tpm table

count = 1
for (abundance in abundances_tsv) {
  tsv_df[count] <- read.csv(abundance, sep = "\t")$tpm
  count <- count + 1
}

transcript_df <- tsv_df

# write.csv(transcript_df, "../transcript_tpm.csv")

## convert transcript id to gene id and aggregate

id_df <- tsv_df
id_df$id <- sub("*\\.[0-9]", "", reduced_id$id)

id_df <- aggregate(id_df[,sapply(id_df,is.numeric)], id_df["id"], sum)
row.names(id_df) <- id_df$id
id_df <- select(id_df, all_of(col))

# write.csv(id_df, "../id_tpm.csv")

## aggregate tpm table to gene-level

tsv_df$symbol <- reduced_gene$symbol
tsv_df <- aggregate(tsv_df[,sapply(tsv_df,is.numeric)], tsv_df["symbol"], sum)
row.names(tsv_df) <- tsv_df$symbol
tsv_df <- select(tsv_df, all_of(col))

gene_df <- tsv_df

# write.csv(gene_df, "../gene_tpm.csv")


## statistical analysis

# hierarchical clustering

# outlier samples?? maybe toss bad outliers
# check with PCA, check number of reads/samples with that sample


gene_df <- as.data.frame(txi_sum$abundance)

dist <- vegdist(t(gene_df), method = "bray")
anova <- anova(betadisper(dist, metadata$Diagnosis))
clus <- hclust(dist, "single",)
plot(clus)

# permanova distances
# deseq finds simple differences, while adonis looks for multivariate differences (differences in genes rely on each other)
# adonis coefficients??
# make captions for all figures, methods for figure generation


permanova <- adonis(t(gene_df) ~ Diagnosis, data = metadata, method = "bray", permutations = 99)
print(as.data.frame(permanova$aov.tab)["Diagnosis", "Pr(>F)"])

coef <- coefficients(permanova)["Diagnosis1",]
top_coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top_coef), horiz = T, las = 1, main = "Most Differential Genes")

bar_data <- as.data.frame(sort(top_coef))
colnames(bar_data) <- "value"
bar_data$labels <- row.names(bar_data)

barplot <- ggplot(data = bar_data, aes(x = reorder(labels, value), y = value, fill = value < 0))

barplot +
  coord_flip() +
  geom_bar(stat="identity") +
  # scale_fill_brewer(palette="Blues") +
  xlab("Gene ID") + 
  ylab("PERMANOVA Linear Model Coefficient") +
  ggtitle("Most Differentially Expressed Genes (Multivariate)") + 
  theme_minimal() +
  theme(legend.position = "none", text = element_text(family = "Helvetica Neue"))

ggsave(filename = "permanova_barplot.png", path = "data/figures/", width = 8, height = 6, device='png', dpi=700)

