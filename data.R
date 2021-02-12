library(Biobase)
library(GEOquery)

protect_soft <- getGEO('GSE109142',GSEMatrix=FALSE)
protect_list <- GSMList(protect_soft)
protect_seq <- list()

count = 1
for (patient in protect_list) {
  protect_seq[[count]] <- read.delim(paste('./data/protect/', list.files('./data/protect')[count], sep=""))
  count <- count + 1
}

risk_soft <- getGEO('GSE117993',GSEMatrix=FALSE)
risk_list <- GSMList(risk_soft)
risk_seq <- list()

count = 1
for (patient in risk_list) {
  risk_seq[[count]] <- read.delim(paste('./data/risk/', list.files('./data/risk')[count], sep=""))
  count <- count + 1
}
