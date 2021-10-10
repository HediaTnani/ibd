library(Biobase)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggdendro)

## read in consistency scores and fluxes

metadata <- read.csv("data/SraRunTable.txt")
# consistency <- read.csv("data/consistency.txt", sep = "\t", header = FALSE)
# consistency <- arrange(consistency, V1, desc())
# rownames(consistency) <- consistency$V1
# consistency <- select(consistency, V2)
# colnames(consistency) <- "consistency"

subsystems <- read.csv("data/subsystems.csv", header = FALSE)

## consistency score dendrogram

# dist <- vegdist(consistency, method = "bray")
# anova <- anova(betadisper(dist, metadata$Diagnosis))
# clus <- ggdendrogram(hclust(dist, "single",))
# clus +
#   ggtitle("GIMME Consistency Scores") +
#   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

## table of fluxes

reactions <- read.table("../Human-GEM-1.9.0/model/reactions1_9.tsv", header=TRUE, sep="\t")
reactions <- select(reactions, rxns, rxnHMR2ID)
colnames(reactions) <- c("mar", "hmr")
flux <- read.table("data/fluxes/SRR6467548.tsv", header = TRUE, sep="\t")
colnames(flux) <- c("hmr", "fluxes")
flux <- select(flux, hmr)
mar <- left_join(flux, reactions, by = "hmr")


fluxes_list <- list.files("data/fluxes")

for (file in fluxes_list){
  if (!exists("fluxes")){
    fluxes <- read.table(paste("data/fluxes/", file, sep = ""), header = TRUE, sep="\t")
    rownames(fluxes) <- mar$hmr
    fluxes <- select(fluxes, "fluxes")
    colnames(fluxes) <- sub('\\.tsv$', '', file)
    fluxes$subsystems <- subsystems$V2
    fluxes$mar <- mar$mar
    fluxes <- fluxes[,c(3,2,1)]
  } else if (exists("fluxes")){
    temp_fluxes <- read.table(paste("data/fluxes/", file, sep = ""), header = TRUE, sep="\t")
    temp_fluxes <- select(temp_fluxes, "fluxes")
    colnames(temp_fluxes) <- sub('\\.tsv$', '', file) 
    fluxes <- cbind(fluxes, temp_fluxes)
    rm(temp_fluxes)
  }
}

# write.csv(fluxes, "data/fluxes.csv")

filtered <- fluxes[rowSums(fluxes[c(3:228)])>0.01 | rowSums(fluxes[c(3:228)]) < -0.01 ,]

write.csv(filtered, "data/filtered_fluxes.csv")

## add subsystems


# permanova <- adonis(t(filtered[c(2:227)]) ~ Diagnosis, data = metadata, method = "bray", permutations = 99)
# print(as.data.frame(permanova$aov.tab)["Diagnosis", "Pr(>F)"])
# 
# coef <- coefficients(permanova)["Diagnosis1",]
# top_coef <- coef[rev(order(abs(coef)))[1:20]]
# par(mar = c(3, 14, 2, 1))
# barplot(sort(top_coef), horiz = T, las = 1, main = "Most Differential Reactions")
