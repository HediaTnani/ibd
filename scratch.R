# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)

# Probiotics intervention example data 
data(peerj32) # Source: https://peerj.com/articles/32/
pseq <- peerj32$phyloseq # Rename the example data

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "group", size = 3)
print(p)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ group,
                    data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])

# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$group))

coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")