### Manhattan plots ###
This script shows how Manhattan plots to represent the SNPs under positive selection
identified by BayeScan, Ohana and PCAdapt were done.

# BayeScan #

setwd("path/to/directory")
library(dplyr)
library(qqman)

# Load data 
bayescan <- read.csv("bayescan_manhattan.csv")

# For SNPs to be considered under diversifying selection in Bayescan,they should meet two criteria:
# qvalue <- 0.05, and ALPHA > 1.0649

# Filter the data to get SNPs under diversifying selection
diversifying_snp <- subset(bayescan, SELECTION == "diversifying")

# Create a vector to store SNPs  under diversifying selection
highlight_snps <- diversifying_snp$SNPs

# Define a Q-value threshold
threshold_qvalue <- 0.05
threshold_logq <- -log10(threshold_qvalue) 

# Add a new column with absolute values of alpha
bayescan$alpha.abs <- abs(bayescan$ALPHA)

# Define alpha threshold
threshold_alpha_min <- 1.06
threshold_alpha_max <- 2.25

# Count the number of SNPs with ALPHA higher than 1.0649
count_alpha_above_threshold <- sum(bayescan$ALPHA > 1.0649, na.rm = TRUE)

# Filter the SNPs with ALPHA higher than 1.0649
snps_above_threshold <- subset(bayescan, ALPHA > 1.0649)

# Manhattan plot using ALPHA values
manhattan(
  bayescan,
  main = "a) BayeScan",
  ylim = c(0, 3),
  chr = "Chr",
  bp = "Position",
  p = "alpha.abs",
  snp = "SNPs",
  ylab = "Alpha",
  col = c("gray20", "gray"),
  chrlabs = NULL,
  suggestiveline = threshold_alpha_min , 
  genomewideline = 5,
  highlight = highlight_snps, 
  logp = FALSE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.1,    #  axis label size
  cex.lab = 1.3,     #  axis title size
  font.lab = 2       #  axis titles bold
)

#################################################################################################################################################
# Ohana #

# Load data 
ohana <- read.csv("manhattan_ohana.csv")

# Filter the data to get SNPs under diversifying selection
diversifying_snp <- subset(ohana, Significance == "Outlier")

# Create vector to store SNPs under diversifying selection
highlight_snps <- diversifying_snp$SNP_ID

# Filter the significant SNPs associated with differences between red and green frogs
significant_redvsgreen <- ohana[ohana$Significance == "Outlier"  & ohana$lle.ratio > 4.9,]

# Create a vector with the SNPs to highlight
highlight_snps <- significant_redvsgreen$SNP_ID

# Define a lle.ratio threshold
threshold_lle.ratio_min <- 4.9
threshold_lle.ratio_max <- 7.5
threshold_loglle.ratio <- log10(threshold_lle.ratio_min) 

# Add the log-transformed lle.ratio to the dataset
ohana$loglle.ratio <- log10(ohana$lle.ratio)

# Count how many SNPs have lle.ratio below the threshold (4.9)
count_snps_below_threshold <- sum(ohana$lle.ratio < threshold_lle.ratio_min, na.rm = TRUE)  ### 95081 SNPs

# Count how many SNPs have lle.ratio above the threshold (4.9)
count_snps_above_threshold <- sum(ohana$lle.ratio > threshold_lle.ratio_min, na.rm = TRUE) ### 910 SNPs

# Count how many SNPs have lle.ratio above the threshold and are marked as Outliers
count_outliers_above_threshold <- sum(ohana$lle.ratio > threshold_lle.ratio_min & ohana$Significance == "Outlier", na.rm = TRUE)
# here the 910 SNP that are significant for differentiating red and green populations

### Manhattan plot using lle.ratio without log10 #####

manhattan(
  ohana,
  main = "b) Ohana ",
  ylim = c(0, 10),
  chr = "Chr",
  bp = "Position",
  p = "lle.ratio",
  snp = "SNP_ID",
  ylab = "Likelihood ratio",
  col = c("gray20", "gray"),
  chrlabs = NULL,
  suggestiveline = threshold_lle.ratio_min , 
  genomewideline = 12,
  highlight = highlight_snps, 
  logp = FALSE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.1,    #  axis label size
  cex.lab = 1.3,     #  axis title size
  font.lab = 2       #  axis titles bold
)

###########################################################################################################################################
# PCAdapt #

# Load data
pcadapt <- read.csv("combined_table_pcadapt.csv")

# Add a new column with absolute values of z.score1
pcadapt$z.score1.abs <- abs(pcadapt$z.score1)

# Filter the significant SNPs associated with PC1
significant_pc1 <- pcadapt[pcadapt$outlier == "Outlier" & !is.na(pcadapt$PC) & pcadapt$PC == 1 & pcadapt$qvalue < 0.05,]
summary(significant_pc1$z.score1.abs)

# Define z.score threshold
threshold_z.score1.abs <- 3.63

# Create a vector with the SNPs to highlight
highlight_snps <- significant_pc1$SNP

# Manhattan plot with absolute values of z-scores, real scale

manhattan(
  pcadapt,   
  main = "c) PCAdapt",
  ylim = c(0, 250),
  chr = "Chr",    # The column with chromosome number
  bp = "POS",     # The column with position
  p = "z.score1.abs",  # The column with absolute Z-scores values
  snp = "SNP",         # The column with SNP names
  ylab = "Absolute Z-scores",  
  col = c("gray20", "gray"),   
  chrlabs = NULL,
  suggestiveline = threshold_z.score1.abs , 
  genomewideline = 300,
  highlight = highlight_snps,
  logp = FALSE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.1,    #  axis label size
  cex.lab = 1.3,     #  axis title size
  font.lab = 2       #  axis titles bold
)
