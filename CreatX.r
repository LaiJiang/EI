
rm(list=ls())

library(xgboost)
library(mlr)
library(pROC)
#######################################################

#load your annottion data raw.annot
load("annot.RData",verbose=TRUE)

test.annot <- raw.annot

test.traits <- unique(test.annot$trait)

snp.dat <- test.annot

#process features
snp.dat$gene.length <- abs(snp.dat$gene.chrom.end - snp.dat$gene.chrom.start) 
snp.dat$snp.tss.dist <- abs(snp.dat$snp.tss.dist)
snp.dat$snp.gene.dist <- abs(snp.dat$snp.gene.dist)
snp.dat$beta <- abs(snp.dat$beta)
snp.dat$z <- abs(snp.dat$z)

###########################################################################
###########################################################################
###########################################################################

source(  "func_summary.r"  )

all.func.names <- c(lsf.str())
all.func.names <- all.func.names[grep("fn.",all.func.names)]  
###########################################################################

snp.collapse <- perpersonsummary(snp.dat,base::mget(all.func.names)  )


Genes.names <- table(snp.dat$gene.name)
Genes.names <- Genes.names[Genes.names>0]
rownames(snp.collapse) <- names(Genes.names)
colnames(snp.collapse) <- all.func.names
###########################################################################

save(snp.collapse , file=paste0(path.dat, test.traits,".RData"))








