### BayeScan analyses ###
This script contains all the steps followed to run BayeScan, from the input file preparation in the adequate format,
to running the software in HPC server, and then analysing results in R. The script for both references was identical, just the input file changes.

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 16.4.2025

help: https://github.com/laurabenestan/Bayescan/blob/master/bayescan.Rmd
help: https://software.bioinformatics.unibe.ch/pgdspider

# Prepare BayeScan input file using sambar R package and PGDSpider #
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.10.txt")
getpackages() ### this is to get all required packages and dependencies sambaR need, and load them

library(adegenet)
library(vcfR)
library (LEA)

# Load data
samples<-read.csv("populations.ingroup.tab", sep="\t",header=T)
importdata(inputprefix="granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK",
           samplefile="populations.ingroup.tab", pop_order=c("DAM","ESQ","SAN","PAL"), 
           colourvector=c("green","springgreen4","coral","red4"),
           geofile="geofile.txt", do_citations=FALSE)

filterdata(indmiss=1,snpmiss=1,min_mac=0,dohefilter=FALSE,min_spacing=0,
           nchroms=NULL,TsTvfilter=NULL)
           
createbayescaninput(allpairwise=FALSE)
# The output file from this step, is then converted to geste format in PGDSpider 3.0.0.0
# Then this file ingeste format will be the input for running BayeScan in Linux server.

# Run BayeScan in Linux server #
#!/bin/bash
#SBATCH --nodes=10  --cpus-per-task=96
#SBATCH --time=48:00:00 
#SBATCH --mem=300G 

# Load necessary modules
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1

threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "path/to/BayeScan software directory"
mkdir output
Bayescan pheno.filter2.bayescan.txt \
-od green_vs_red \
-threads $SLURM_CPUS_PER_TASK -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10

# Analize BayeScan output file in R #

setwd("path/to/directory")
library(ggplot2)

# Load data
bayescan=read.table("pheno.filter2.bayescan_fst.txt") 
# The first column of the `bayescan` object is a SNP ID. 
#The next three (prob, log10(P0), and qval) are related to the test of local adaptation considering the logarithm of the posterior odds - log10(PO)
# and the q-value for the model with selection. The fifth column gives the size of the locus-specific effect (alpha parameter). 
# The last one provides the locus-specific FST averaged over all populations.

SNP=read.table("pheno.filter2.letter.map",header=FALSE)
# Extract SNP_IDs in the 2nd column for the SNP_ID table
library(dplyr) 
SNP %>% select(1,2,4)
# Creating SNP_ID by joining col 1 (chromosome) and col4 (position)
SNP_ID <- paste(SNP$V1, SNP$V4, sep="_") 

# Merge the name of the outliers with the results from BayeScan. 
bayescan=cbind(SNP$V2, bayescan) 
bayescan=cbind(SNP_ID, bayescan) 
colnames(bayescan)=c("SNP_ID","SNPs","PROB","LOG_PO","Q_VALUE","ALPHA","FST")

# Change the value of the Q_VALUE column: 0 == 0.0001. 
attach(bayescan)
class(bayescan$Q_VALUE) 
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

# Round the values 
bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST,6))

# Add a column for the type of selection grouping based on a Q-VALUE < 0.05
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

# Save the results of the SNPs potentially under positive and balancing selection (qvalue < 0.05)
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

# Check the number of SNPs belonging to each category
xtabs(data=bayescan, ~SELECTION) 

# Write the results of the SNPs potentially under selection (qvalue < 0.05). 
write.table(neutral, "granulifera_neutral_SNPs.txt", row.names=F, quote=F) 
write.table(balancing, "granulifera_balancing_SNPs.txt", row.names=F, quote=F) 
write.table(positive, "granulifera_positive_SNPs.txt", row.names=F, quote=F)

# Transformation Log of the Q-value in order to create ggplot graph
range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE)

# Create title for the ggplot graph 
x_title="Log(q-value)" 
y_title="Fst"

# Make the ggplot graph 
graph_1<-ggplot(bayescan,aes(x=LOG10_Q,y=FST, label=bayescan$POS)) 
graph_1+geom_point(aes(color=bayescan$SELECTION), pch=20, size=2)+ 
    #geom_text()+ 
    scale_color_manual(name="Selection",values=c("black","red", "gray"))+ 
    labs(x=x_title)+ 
    labs(y=y_title)+ 
    ggtitle("Bayescan plot Oophaga granulifera")+
    theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+ 
    theme(axis.text.x=element_text(colour="black"))+ 
    theme(axis.text.y=element_text(colour="black",size=12))+ 
    theme(axis.text.x=element_text(colour="black",size=12))+ 
    theme(panel.border = element_rect(colour="black", fill=NA, size=3),  
          axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold")) +
    theme_classic()
  
# Save the figure
ggsave("bayescan_granulifera_genome.pdf", dpi=600, width=5, height=5) 
ggsave("bayescan_granulifera_genome.jpeg", dpi=600, width=7, height=5) 
dev.off()




 
