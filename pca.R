library(Biobase)
library(dplyr)
library(ggplot2)
library(vegan)
library(randomForest)

pca_out <- read.csv("data/gimme_pca.csv")

rownames(pca_out) <- pca_out$X
pca_out <- pca_out[,2:4]

# sort pca by groups

count <- 1
for (i in pca_out$PC1) {
  if (i < 0) {
    pca_out$group[count] <- 1
  } else if (i >= 0 & pca_out$PC2[count] > -250) {
    pca_out$group[count] <- 2
  } else if (i >= 0 & pca_out$PC2[count] < -250) {
    pca_out$group[count] <- 3
  } else {
    next
  }
  count <- count + 1
}

## folate pca violin plot


fluxes_old$diagnosis <- as.factor(metadata$Diagnosis)
fluxes_old <- fluxes_old %>%
  select(diagnosis, everything())


fol <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
fol$subsystem <- description
fol$mar <- feature_list$mar
fol <- filter(fol, subsystem == "Folate metabolism")
rownames(fol) <- fol$mar
fol <- fol[, 1:(ncol(fol)-2)]
fol <- as.data.frame(t(fol))
fol$diagnosis <- as.character(pca_out$group)
fol <- select(fol, c(-"MAR08144",-"MAR03929",-"MAR04333",-"MAR04655",-"MAR04654",-"MAR04332"))
fol_remove <-c("SRR6467726","SRR6467645","SRR6467687","SRR6467734")
fol <- fol[!(row.names(fol) %in% fol_remove),]
fol_melt <- fol %>% 
  tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
  mutate(text = reorder(text, value, mean))

# v_test <- ggplot(fol, aes(x=MAR03929, y=diagnosis, fill=diagnosis))
# v_test +
#   geom_violin()

v1 <- ggplot(fol_melt, aes(x=text, y=value, fill=diagnosis))

v1 +
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  ylab("Flux")+
  xlab("Reactions")+
  ggtitle("Folate Metabolism Fluxes") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  coord_flip() + 
  scale_fill_discrete(name = "PCA Groups") +
  theme(legend.position = c(0.85, 0.2)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank())

ggsave(filename = "violin_fol.png", path = "data/figures/", width = 6, height = 8, device='png', dpi=700)

## pentose cycle violin plot

pentose <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
pentose$subsystem <- description
pentose$mar <- feature_list$mar
# pentose <- filter(pentose, subsystem == "Transport reactions")
pentose <- filter(pentose, mar == "MAR00788")
rownames(pentose) <- pentose$mar
pentose <- pentose[, 1:(ncol(pentose)-2)]
pentose <- as.data.frame(t(pentose))
pentose$diagnosis <- as.character(pca_out$group)
pentose_c <- filter(pentose, diagnosis == "Control")
pentose_u <- filter(pentose, diagnosis == "Ulcerative Colitis")
# pentose <- select(pentose, c(-"MAR08653",-"MAR04304",-"MAR04567"))
pentose <- select(pentose, c("MAR00788","diagnosis"))
pentose_melt <- pentose %>% 
  tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
  mutate(text = reorder(text, value, mean))

# v4 <- ggplot(pentose_melt, aes(x=text, y=value, fill=diagnosis))
v4 <- ggplot(pentose, aes(x=MAR00788, y=diagnosis, fill=diagnosis))

v4 +
  geom_violin(show.legend = FALSE)+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA,show.legend = FALSE)+
  # ylab("Flux")+
  ylab("MAR00788")+
  # xlab("Reactions")
  xlab("Flux")+
  # ggtitle("Pentose Phosphate Pathway Reaction Fluxes") +
  ggtitle("Acyl-CoA Dehydrogenase Flux (MAR00788)") +
  theme_minimal()+
  scale_y_discrete(labels = c("Group 1", "Group 2", "Group 3"))+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9),show.legend = FALSE) + 
  # coord_flip() + 
  # scale_fill_discrete(name = "PCA Group") +
  # theme(legend.position = c(0.85, 0.4)) +
  # theme(legend.position = c(0.2, 0.85)) +
  # theme(legend.background = element_rect(fill="white",
  #                                        size=0.2, linetype="solid", 
  #                                        colour ="black")) +
  theme(plot.title = element_text(hjust = 1.1), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank(), axis.title.y=element_blank())

ggsave(filename = "violin_sh.png", path = "data/figures/", width = 4, height = 2, device='png', dpi=1400)


# ml

metadata <- read.csv("data/SraRunTable.txt")

fluxes <- read.csv("data/filtered_fluxes.csv")
rownames(fluxes) <- fluxes$X
reaction <- fluxes$X
description <- fluxes$subsystems
fluxes <- fluxes[,4:ncol(fluxes)]
fluxes <- as.data.frame(t(fluxes))
fluxes[is.na(fluxes)] <- 0
fluxes$diagnosis <- as.factor(metadata$Diagnosis)

fluxes$run <- metadata$Run
fluxes <- filter(fluxes, run != "SRR6467726")
new_d <- fluxes$diagnosis
fluxes <- select(fluxes, c(-"run",-"diagnosis"))

fluxes_old <- fluxes
colnames(fluxes) <- sub('\\[e]$', '', colnames(fluxes))
colnames(fluxes) <- paste('reaction_', colnames(fluxes), sep='')

fluxes$group <- as.factor(pca_out$group)
fluxes_old$group <- as.factor(pca_out$group)

fluxes <- fluxes %>%
  select(group, everything())

fluxes_old <- fluxes_old %>%
  select(group, everything())



# fluxes_old <- as.data.frame(t(fluxes_old))

training <- fluxes[c(1:15, 21:176),]
validation <- fluxes[c(16:20, 177:226),]

forest <- randomForest(group ~ ., data = fluxes, ntree=400, mtry=25, importance=TRUE)

print(forest)

# roc_prediction <- as.data.frame(predict(forest, validation[,-650]))

feature_importance <- as.data.frame(importance(forest))
feature_importance$description <- description
feature_importance$hmr <- reaction

reactions <- read.table("../Human-GEM-1.9.0/model/reactions1_9.tsv", header=TRUE, sep="\t")
hmr <- select(reactions, rxns, rxnHMR2ID)
colnames(hmr) <- c("mar", "hmr")
recon <- select(reactions, rxns, rxnRecon3DID)
colnames(recon) <- c("mar", "hmr")
features_mar <- select(feature_importance, hmr)
features_mar <- left_join(features_mar, hmr, by = "hmr")
features_recon <- features_mar %>%
  filter_all(any_vars(is.na(.)))
features_mar <- na.omit(features_mar)
features_recon <- select(features_recon, hmr)
features_recon <- left_join(features_recon, recon, by = "hmr")
features <- rbind(features_mar, features_recon)
feature_importance <- left_join(feature_importance, features, by = "hmr")


fluxes_1 <- fluxes_old[fluxes_old$group == 1, 2:ncol(fluxes_old)]
means <- as.data.frame(colMeans(fluxes_1))
fluxes_2 <- fluxes_old[fluxes_old$group == 2, 2:ncol(fluxes_old)]
means$means_2 <- colMeans(fluxes_2)
fluxes_3 <- fluxes_old[fluxes_old$group == 3, 2:ncol(fluxes_old)]
means$means3 <-colMeans(fluxes_3)
means$hmr <- feature_importance$hmr
means$mar <- feature_importance$mar
colnames(means) <- c("means_1", "means_2", "means_3", "hmr", "mar")

# write.csv(means, "data/feature_means_group.csv")

# rownames(feature_importance) <- reaction

# write.csv(feature_importance, "data/feature_importance_group.csv")



# permanova

# data <- fluxes_old[,2:ncol(fluxes_old)]
# data <- mutate_all(fluxes_old, function(x) as.numeric(as.character(x)))
# # data$group <- as.factor(pca_out$group)
# 
# permanova <- adonis(data~group, data = pca_out, method = "euclidean", permutations = 99)
# print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
# 
# coef <- coefficients(permanova)["group1",]
# p_results <- as.data.frame(reaction)
# p_results$description <- as.data.frame(description)
# p_results$coef <- as.data.frame(coef)
# top_coef <- coef[rev(order(abs(coef)))[1:20]]
# par(mar = c(3, 14, 2, 1))
# barplot(sort(top_coef), horiz = T, las = 1, main = "Most Differential Reactions")
# 
# bar_data <- as.data.frame(sort(top_coef))
# colnames(bar_data) <- "value"
# bar_data$labels <- row.names(bar_data)
# 
# barplot <- ggplot(data = bar_data, aes(x = reorder(labels, value), y = value, fill = value < 0))
# 
# barplot +
#   coord_flip() +
#   geom_bar(stat="identity") +
#   # scale_fill_brewer(palette="Blues") +
#   xlab("Reaction ID") +
#   ylab("PERMANOVA Linear Model Coefficient") +
#   ggtitle("Most Differentially Regulated Reactions") +
#   theme_minimal() +
#   theme(legend.position = "none", text = element_text(family = "Helvetica Neue"))
# 
# ggsave(filename = "gimme_barplot.png", path = "data/figures/", width = 8, height = 6, device='png', dpi=700)
# 
# ## wald test
# 
# lm <- lm(group ~ ., data = fluxes_old)
# 
# wald <- wald.test(lm$model, lm$coefficients)
# 


