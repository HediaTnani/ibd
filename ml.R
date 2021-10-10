library(Biobase)
library(dplyr)
library(party)
library(randomForest)
require(caTools)
library(ROCR)
library(ggplot2)
library(vegan)
library(devtools)
library(ggbiplot)

metadata <- read.csv("data/SraRunTable.txt")

fluxes <- read.csv("data/filtered_fluxes.csv")
rownames(fluxes) <- fluxes$X
reaction <- fluxes$X
description <- fluxes$subsystems
fluxes <- fluxes[,4:ncol(fluxes)]
fluxes <- as.data.frame(t(fluxes))
fluxes$diagnosis <- as.factor(metadata$Diagnosis)

fluxes$run <- metadata$Run
# fluxes <- filter(fluxes, run != "SRR6467726")
new_d <- fluxes$diagnosis
fluxes <- select(fluxes, c(-"run",-"diagnosis"))

fluxes_old <- fluxes
colnames(fluxes) <- sub('\\[e]$', '', colnames(fluxes))
colnames(fluxes) <- paste('reaction_', colnames(fluxes), sep='')
fluxes[is.na(fluxes)] <- 0

fluxes$diagnosis <- as.factor(new_d)

training <- fluxes[c(1:15, 21:176),]
validation <- fluxes[c(16:20, 177:226),]

## remission score

remission_na <- fluxes_old
remission_na$score <- as.factor(metadata$week_4_remission)

remission_na <- remission_na %>%
  select(score, everything())

remission <- na.omit(remission_na)

fluxes_yes <- remission[remission$score == "Yes", 2:ncol(remission)]
means_rem <- as.data.frame(colMeans(fluxes_yes))
fluxes_no <- remission[remission$score == "No", 2:ncol(remission)]
means_rem$means_1 <- colMeans(fluxes_no)

colnames(means_rem) <- c("means_yes", "means_no")

feature_list <- as.data.frame(as.character(reaction))
colnames(feature_list) <- "hmr"
reactions <- read.table("../Human-GEM-1.9.0/model/reactions1_9.tsv", header=TRUE, sep="\t")
hmr <- select(reactions, rxns, rxnHMR2ID)
colnames(hmr) <- c("mar", "hmr")
recon <- select(reactions, rxns, rxnRecon3DID)
colnames(recon) <- c("mar", "hmr")
features_mar <- left_join(feature_list, hmr, by = "hmr")
features_recon <- features_mar %>%
  filter_all(any_vars(is.na(.)))
features_mar <- na.omit(features_mar)
features_recon <- select(features_recon, hmr)
features_recon <- left_join(features_recon, recon, by = "hmr")
features <- rbind(features_mar, features_recon)
feature_list <- left_join(feature_list, features, by = "hmr")

means_rem$subsystem <- description
means_rem$hmr <- feature_list$hmr
means_rem$mar <- feature_list$mar

## histology score

hist_fluxes_na <- fluxes_old
hist_fluxes_na$score <- as.factor(metadata$histology_severity_score)

hist_fluxes_na <- hist_fluxes_na %>%
  select(score, everything())

hist_fluxes <- na.omit(hist_fluxes_na)

# hist_fluxes <- as.data.frame(t(hist_fluxes))

fluxes_0 <- hist_fluxes[hist_fluxes$score == 0, 2:ncol(hist_fluxes)]
means_hist <- as.data.frame(colMeans(fluxes_0))
fluxes_1 <- hist_fluxes[hist_fluxes$score == 1, 2:ncol(hist_fluxes)]
means_hist$means_1 <- colMeans(fluxes_1)
fluxes_2 <- hist_fluxes[hist_fluxes$score == 2, 2:ncol(hist_fluxes)]
means_hist$means_2 <- colMeans(fluxes_2)
fluxes_3 <- hist_fluxes[hist_fluxes$score == 3, 2:ncol(hist_fluxes)]
means_hist$means3 <-colMeans(fluxes_3)
fluxes_4 <- hist_fluxes[hist_fluxes$score == 4, 2:ncol(hist_fluxes)]
means_hist$means4 <-colMeans(fluxes_4)
# means$hmr <- feature_importance$hmr
# means$mar <- feature_importance$mar
colnames(means_hist) <- c("means_0", "means_1", "means_2", "means_3", "means_4")

means_hist$subsystem <- description
means_hist$hmr <- feature_list$hmr
means_hist$mar <- feature_list$mar

## pucai score

pucai_fluxes_na <- fluxes_old
pucai_fluxes_na$score <- as.factor(pucai$group)

pucai_fluxes_na <- pucai_fluxes_na %>%
  select(score, everything())

pucai_fluxes <- na.omit(pucai_fluxes_na)

# pucai_fluxes <- as.data.frame(t(hist_fluxes))

fluxes_0 <- pucai_fluxes[pucai_fluxes$score == 0, 2:ncol(pucai_fluxes)]
means_pucai <- as.data.frame(colMeans(fluxes_0))
fluxes_1 <- pucai_fluxes[pucai_fluxes$score == 1, 2:ncol(pucai_fluxes)]
means_pucai$means_1 <- colMeans(fluxes_1)
fluxes_2 <- pucai_fluxes[pucai_fluxes$score == 2, 2:ncol(pucai_fluxes)]
means_pucai$means_2 <- colMeans(fluxes_2)
fluxes_3 <- pucai_fluxes[pucai_fluxes$score == 3, 2:ncol(pucai_fluxes)]
means_pucai$means3 <-colMeans(fluxes_3)
fluxes_4 <- pucai_fluxes[pucai_fluxes$score == 4, 2:ncol(pucai_fluxes)]
means_pucai$means4 <-colMeans(fluxes_4)

colnames(means_pucai) <- c("means_0", "means_1", "means_2", "means_3", "means_4")

means_pucai$subsystem <- description
means_pucai$hmr <- feature_list$hmr
means_pucai$mar <- feature_list$mar

write.csv(means_pucai, "data/pucai_means.csv")

## enrichment analysis

# fishers <- fisher.test(fluxes_old)

## pca

# do enrichment analysis to identify subsystems in 3 groups
# each group may be led by certain significant pathways
# reactions and subsystems -> ratio of categories in patients
# fisher's exact?

pucai <- as.data.frame(metadata$PUCAI)
pucai[is.na(pucai)] <- 0
colnames(pucai) <- "PUCAI"

count <- 1
for (i in pucai$PUCAI) {
  if (i == 0) {
    pucai$group[count] <- 0
  } else if (i >= 1 & i < 22) {
    pucai$group[count] <- 1
  } else if (i >= 22 & i < 43) {
    pucai$group[count] <- 2
  } else if (i >= 43 & i < 64) {
    pucai$group[count] <- 3
  } else if (i >= 64) {
    pucai$group[count] <- 4
  } 
  count <- count + 1
}

no_transport <- as.data.frame(t(fluxes_old))
no_transport$subsystem <- description
no_transport <- filter(no_transport, subsystem != "Transport reactions" & subsystem != "Exchange/demand reactions")
desc_nt <- no_transport$subsystem
no_transport <- select(no_transport, -"subsystem")
no_transport <- as.data.frame(t(no_transport))

pca <- prcomp(fluxes_old)
# pca <- prcomp(no_transport)
# summary(pca)

# plot(pca$x[,1], pca$x[,2])
pca_out <- as.data.frame(pca$x)
pca_out$diagnosis <- as.character(new_d)
# pca_out$diagnosis <- as.character(metadata$histology_severity_score)
# pca_out$diagnosis <- as.factor(metadata$week_4_remission)
# pca_out$diagnosis <- as.factor(pucai$group)
pca_out <- pca_out %>%
  select(diagnosis, everything())

write.csv(pca_out, "data/gimme_pca.csv")

p <- ggplot(pca_out, aes(x=PC1, y=PC2, color=diagnosis))
p +
  geom_point(aes(colour = diagnosis), alpha = 0.5,subset(pca_out,diagnosis != "Control")) +
  geom_point(aes(colour = diagnosis), alpha = 0.5,subset(pca_out,diagnosis != "Ulcerative Colitis")) +
  ggtitle("Patient Fluxes") +
  theme_minimal() +
  xlab("PC1 (32.90%)") + 
  ylab("PC2 (22.87%)") + 
  labs(colour = "Diagnosis")  +
  # labs(colour = "PUCAI Group")  +
  theme(legend.position = c(0.225, 0.225)) +
  # theme(legend.position = c(0.225, 0.25)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))
  

ggsave(filename = "pca.png", path = "data/figures/", width = 6, height = 4.5, device='png', dpi=1400)

## random forest classifier

# fluxes$pucai <- as.factor(pucai$group)


colnames(no_transport) <- sub('\\[e]$', '', colnames(no_transport))
colnames(no_transport) <- paste('reaction_', colnames(no_transport), sep='')
no_transport[is.na(no_transport)] <- 0

no_transport$diagnosis <- as.factor(metadata$Diagnosis)

# forest <- randomForest(pucai ~ ., data = fluxes, ntree=300, mtry=25, importance=TRUE)
forest <- randomForest(diagnosis ~ ., data = no_transport, ntree=300, mtry=25, importance=TRUE, classwt = c(0.9, 0.1))

# class weight parameter (unequal group sizes) penalizes picking larger group

print(forest)

## feature importance

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

# write.csv(feature_importance, "data/feature_importance.csv")

## flux means

fluxes_old$diagnosis <- as.factor(metadata$Diagnosis)
fluxes_old <- fluxes_old %>%
  select(diagnosis, everything())

fluxes_control <- fluxes_old[fluxes_old$diagnosis == "Control", 2:ncol(fluxes_old)]
means <- as.data.frame(colMeans(fluxes_control))
fluxes_uc <- fluxes_old[fluxes_old$diagnosis == "Ulcerative Colitis", 2:ncol(fluxes_old)]
means$means_2 <- colMeans(fluxes_uc)

control_means <- means[,1]
ibd_means <-  means[,2]
fold_change <- ibd_means-control_means
difference <- scale(fold_change)
means$fold_change <- difference

means$hmr <- feature_list$hmr
means$mar <- feature_list$mar
colnames(means) <- c("means_control", "means_uc", "fold_change","hmr", "mar")

# write.csv(means, "data/feature_means.csv")

##glycolysis violin plot

library(grid)
library(gridExtra)
library(tidyverse)

glycolysis <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
glycolysis$subsystem <- description
glycolysis$mar <- feature_list$mar
glycolysis <- filter(glycolysis, subsystem == "Glycolysis / Gluconeogenesis")
rownames(glycolysis) <- glycolysis$mar
glycolysis <- glycolysis[, 1:(ncol(glycolysis)-2)]
glycolysis <- as.data.frame(t(glycolysis))
glycolysis$diagnosis <- metadata$Diagnosis
glycolysis_c <- filter(glycolysis, diagnosis == "Control")
glycolysis_u <- filter(glycolysis, diagnosis == "Ulcerative Colitis")
glycolysis <- select(glycolysis, c(-"MAR04355",-"MAR04358",-"MAR04371"))
glycolysis_melt <- glycolysis %>% 
  tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
  mutate(text = reorder(text, value, mean))

v1 <- ggplot(glycolysis_melt, aes(x=text, y=value, fill=diagnosis))

v1 +
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  ylab("Flux")+
  xlab("Reactions")+
  ggtitle("Glycolysis Reaction Fluxes") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  scale_x_discrete(labels = c("MAR04373\n(Glyceraldehyde-\n3-phosphate\nDehydrogenase)","MAR04375\n(Aldolase)","MAR04365\n(Phosphoglycerate\nMutase)","MAR04396\n(Phosphoglucomutase)","MAR04391\n(Triosephosphate\nIsomerase)","MAR04394\n(Hexokinase)","MAR04381\n(Glucose-6-phosphate\nIsomerase)","MAR04363\n(Enolase)","MAR04369\n(Phosphoglycerate\nKinase)")) +
  coord_flip() + 
  scale_fill_discrete(name = "Diagnosis") +
  theme(legend.position = c(0.775, 0.175)) +
  theme(legend.background = element_rect(fill="white",
                                             size=0.2, linetype="solid", 
                                             colour ="black")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank())

ggsave(filename = "violin_glycolysis1.png", path = "data/figures/", width = 6, height = 8, device='png', dpi=700)

## oxidative phosphorylation violin plot 

oxi <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
oxi$subsystem <- description
oxi$mar <- feature_list$mar
oxi <- filter(oxi, subsystem == "Oxidative phosphorylation")
rownames(oxi) <- oxi$mar
oxi <- oxi[, 1:(ncol(oxi)-2)]
oxi <- as.data.frame(t(oxi))
oxi$diagnosis <- metadata$Diagnosis
oxi_c <- filter(oxi, diagnosis == "Control")
oxi_u <- filter(oxi, diagnosis == "Ulcerative Colitis")
oxi <- select(oxi, c(-"MAR06914"))
oxi_melt <- oxi %>% 
  tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
  mutate(text = reorder(text, value, mean))

v2 <- ggplot(oxi_melt, aes(x=text, y=value, fill=diagnosis))

v2 +
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  ylab("Flux")+
  xlab("Reactions")+
  ggtitle("Oxidative Phosphorylation Fluxes") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  coord_flip() + 
  scale_fill_discrete(name = "Diagnosis") +
  theme(legend.position = c(0.85, 0.2)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank())

ggsave(filename = "violin_oxi.png", path = "data/figures/", width = 6, height = 8, device='png', dpi=700)

## tca cycle violin plot

tca <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
tca$subsystem <- description
tca$mar <- feature_list$mar
tca <- filter(tca, subsystem == "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")
rownames(tca) <- tca$mar
tca <- tca[, 1:(ncol(tca)-2)]
tca <- as.data.frame(t(tca))
tca$diagnosis <- metadata$Diagnosis
tca_c <- filter(tca, diagnosis == "Control")
tca_u <- filter(tca, diagnosis == "Ulcerative Colitis")
tca <- select(tca, c(-"MAR01453",-"MAR00710",-"MAR04408",-"MAR04589",-"MAR04147"))
tca_melt <- tca %>% 
  tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
  mutate(text = reorder(text, value, mean))

v3 <- ggplot(tca_melt, aes(x=text, y=value, fill=diagnosis))

v3 +
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  ylab("Flux")+
  xlab("Reactions")+
  ggtitle("TCA Cycle Fluxes") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  coord_flip() + 
  scale_fill_discrete(name = "Diagnosis") +
  theme(legend.position = c(0.2, 0.85)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank())

ggsave(filename = "violin_tca.png", path = "data/figures/", width = 6, height = 8, device='png', dpi=700)

## pentose cycle violin plot

pentose <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
pentose$subsystem <- description
pentose$mar <- feature_list$mar
# pentose <- filter(pentose, subsystem == "Transport reactions")
pentose <- filter(pentose, mar == "MAR00788")
rownames(pentose) <- pentose$mar
pentose <- pentose[, 1:(ncol(pentose)-2)]
pentose <- as.data.frame(t(pentose))
pentose$diagnosis <- metadata$Diagnosis
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
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  # ylab("Flux")+
  ylab("MAR00788")+
  # xlab("Reactions")
  xlab("Flux")+
  # ggtitle("Pentose Phosphate Pathway Reaction Fluxes") +
  ggtitle("Acyl-CoA Dehydrogenase Flux (MAR00788)") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  # coord_flip() + 
  scale_fill_discrete(name = "Diagnosis") +
  # theme(legend.position = c(0.85, 0.4)) +
  theme(legend.position = "none") +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  theme(plot.title = element_text(hjust = 1.1), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank(), axis.title.y=element_blank())

ggsave(filename = "violin_acd.png", path = "data/figures/", width = 4, height = 2, device='png', dpi=1400)

## fatty cycle violin plot

fatty <- as.data.frame(t(fluxes_old[,2:ncol(fluxes_old)]))
fatty$subsystem <- description
fatty$mar <- feature_list$mar
# fatty <- filter(fatty, subsystem == "Tryptophan metabolism")
fatty <- filter(fatty, mar == "MAR00788")
rownames(fatty) <- fatty$mar
fatty <- fatty[, 1:(ncol(fatty)-2)]
fatty <- as.data.frame(t(fatty))
fatty$diagnosis <- as.factor(pucai$group)
# fatty_c <- filter(fatty, diagnosis == "Control")
# fatty_u <- filter(fatty, diagnosis == "Ulcerative Colitis")
# fatty <- select(fatty, c(-"MAR00788"))
# fatty_melt <- fatty %>% 
#   tidyr::pivot_longer(-diagnosis, names_to = "text") %>% 
#   mutate(text = reorder(text, value, mean))

v5 <- ggplot(fatty, aes(x=MAR00788, y=diagnosis, fill=diagnosis))

v5 +
  geom_violin()+
  geom_boxplot(width=0.05, position=position_dodge(width=.9),outlier.shape = NA)+
  ylab("Flux")+
  xlab("Reactions")+
  ggtitle("Acyl-CoA Dehydrogenase Flux (MAR00788)") +
  theme_minimal()+
  stat_summary(fun=mean, geom="point", size=1, color="red", position=position_dodge(width=.9)) + 
  # coord_flip() + 
  scale_fill_discrete(name = "PUCAI",labels = c("Control", "0-21","22-42","43-63","64-85")) +
  theme(legend.position = c(0.8, 0.275)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.2, linetype="solid", 
                                         colour ="black")) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))+
  theme(panel.grid.major.y = element_blank(), axis.title.y=element_blank())

ggsave(filename = "violin_acd.png", path = "data/figures/", width = 4, height = 3, device='png', dpi=1400)

## volcano plot

volcano_data <- as.data.frame(difference)
colnames(volcano_data) <- "difference"
volcano_data$importance <- feature_importance$MeanDecreaseGini
volcano_data$reaction <- reaction
volcano_data$description <- description
rownames(volcano_data) <- reaction
volcano_data$label[1] <- ""

count <- 1
for (i in volcano_data$difference) {
  if (i > 2.3 || i < -2.3) {
    if (volcano_data$importance[count] > 0.5) {
      print(count)
      volcano_data$label[count] <- volcano_data$reaction[count]
    }
  } 
  count <- count + 1
}

volcano <- ggplot(volcano_data, aes(x=difference, y=feature_importance$MeanDecreaseGini))
volcano + 
  geom_text(aes(label = label), hjust = 0, vjust = 0) +
  geom_point(aes(alpha = 0.5), show.legend = FALSE) +
  ggtitle("RF Feature Importance and Flux Difference") + 
  theme_minimal() +
  xlab("Flux Difference Between Groups") + 
  ylab("Gini Index") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))

ggsave(filename = "feature_importance.png", path = "data/figures/", width = 8, height = 6, device='png', dpi=700)

## permanova

permanova <- adonis((fluxes_old[,1:649]) ~ Diagnosis, data = metadata, method = "euclidean", permutations = 99)
print(as.data.frame(permanova$aov.tab)["Diagnosis", "Pr(>F)"])

coef <- coefficients(permanova)["Diagnosis1",]
p_results <- as.data.frame(reaction)
p_results$description <- as.data.frame(description)
p_results$coef <- as.data.frame(coef)
top_coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top_coef), horiz = T, las = 1, main = "Most Differential Reactions")

p_out <- as.data.frame(coef)
p_out$description <- description
write.csv(p_out, "data/gimme_permanova.csv")

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
  
## roc curve

roc_prediction <- predict(forest, validation[,-650], type="prob")
classes <- levels(validation$diagnosis)

pretty_colours <- c("#F8766D","#00BA38","#619CFF")

for (i in 1:2)
{
  # Define which observations belong to class[i]
  true_values <- ifelse(validation[,650]==classes[i],1,0)
  # Assess the performance of classifier for class[i]
  pred <- prediction(roc_prediction[,i],true_values)
  perf <- performance(pred, "tpr", "fpr")
  if (i==1)
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i]) 
  }
  else
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i],add=TRUE) 
  }
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  # print(auc.perf@y.values)
}
