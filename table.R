library(Biobase)
library(annotables)
library(dplyr)
library(vegan)
library(data.table)

## produce transcript-level tpm table

ttg <- read.table("data/tx2gene.txt", sep = "\t", col.names = c("target_id", "ens_gene", "ext_gene"))

metadata <- read.csv("../protect/SraRunTable.txt")

tsv_folder <- "../protect/tsv"
tsv_subfolders <- list.dirs(tsv_folder, recursive = TRUE)[-1]
abundances_tsv <- file.path(tsv_subfolders, "abundance.tsv")
abundances_h5 <- file.path(tsv_subfolders, "abundance.h5")

row <- read.csv(abundances_tsv[[1]], sep = "\t")$target_id
col <- metadata$Run
df_row <- data.frame(row)
ttg_id <- select(ttg, c("target_id", "ext_gene"))

names(ttg_id) <- c("row", "symbol")
ttg_id$symbol[[187688]] <- "5-nucleotidase"
ttg_id$symbol[[187461]] <- "Wilms"
reduced <- merge(df_row, ttg_id, all.x=TRUE)  

tsv_df <- data.frame(matrix(ncol=226, nrow=199240, dimnames=list(row, col)))

count = 1
for (abundance in abundances_tsv) {
  tsv_df[count] <- read.csv(abundance, sep = "\t")$tpm
  count <- count + 1
}

transcript_df <- tsv_df

# write.csv(transcript_df, "../transcript_tpm.csv")

## aggregate tpm table to gene-level

tsv_df$symbol <- reduced$symbol
tsv_df <- aggregate(tsv_df[,sapply(tsv_df,is.numeric)], tsv_df["symbol"], sum)
row.names(tsv_df) <- tsv_df$symbol
tsv_df <- select(tsv_df, all_of(col))

gene_df <- tsv_df

# write.csv(gene_df, "../gene_tpm.csv")

## statistical testing

# anova and hierarchical clustering

dist <- vegdist(t(gene_df), method = "euclidean")
anova <- anova(betadisper(dist, metadata$Diagnosis))
clus <- hclust(dist, "single",)
plot(clus)

# transcript permanova

permanova <- adonis(t(transcript_df) ~ Diagnosis, data = metadata, method = "bray", permutations = 99)
print(as.data.frame(permanova$aov.tab)["Diagnosis", "Pr(>F)"])

coef <- coefficients(permanova)["Diagnosis1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top Transcripts")

# gene permanova

permanova <- adonis(t(gene_df) ~ Diagnosis, data = metadata, method = "euclidean", permutations = 99)
print(as.data.frame(permanova$aov.tab)["Diagnosis", "Pr(>F)"])

coef <- coefficients(permanova)["Diagnosis1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top Genes")



