### Selection scan with PCAdapt ### 
The script for both references is identical, just the input file changes.

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 16.4.2025

help:  https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

setwd("path/to/directory")

library(vcfR)
library(pcadapt)
library(ggplot2)
library(qvalue)

# Load data
bed<-read.pcadapt("granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK.bed",
                  type = "bed",type.out = "bed")
pops <- read.delim("populations.ingroup.tab",header=TRUE, stringsAsFactors=FALSE)
population_colors <- c("darkgreen", "lightgreen", "red", "darkred")

x <- pcadapt(input = bed, K = 10,ploidy=2) # k is the number of principal components to be retained, initially I put 10 to explore the data
plot(x, option = "screeplot", K = 10)
plot(x, option = "screeplot", K = 5)
plot(x, option = "scores", pop = pops$pop)
plot(x, option = "scores", i = 3, j = 4, pop = pops$pop)
plot(x, option = "scores", i = 5, j = 6, pop = pops$pop)

# From the screeplot and the scores plots of various PCs it seems that retaining 3 PCs is the best option
x <- pcadapt(input = bed, K = 3)
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")   
plot(x, option = "qqplot")
plot(x, option = "manhattan")
summary(x)
removed_snps <- which(is.na(x$stat))

# Choosing a cutoff for outlier detection 
# Transform p-values into q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers) 
write.csv(outliers, "PCAdapt_outliers.csv")

# Associate PCs and outliers
snp_pc <- get.pc(x, outliers)  
write.csv(snp_pc, "PCAdapt.snps-pc.csv")

# Associate outliers to all statistics provided by PCAdapt
snp_data<-read.delim("granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK.bim",header=FALSE, stringsAsFactors=FALSE)
x$qvalue<- qval
results<-cbind(x$loadings, x$zscores,x$af,x$maf,x$chi2.stat,x$stat,x$pvalues,x$qvalue,snp_data[1:6])
names(results)<-c("PC1","PC2","PC3", "z.score1","z.score2",
                  "z.score3","af", "maf","chi2.stat","stat","pvalues", 
                  "qvalue", "CHROM","SNP","MISS","POS","ref","alt")
write.csv(results,"ALL_SNPs_statistics.csv")

# Extract significant snps
signif<-results[results$qvalue<alpha,]
signif$assoc.PC<-as.numeric(snp_pc$PC)
length(signif$qval)
write.csv(signif,"PCAdapt_significant_SNPS.csv")
