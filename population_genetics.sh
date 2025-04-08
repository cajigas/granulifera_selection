Population genetics analyses were performed using the R package sambar and included:
1- Population structure
2- Population genetic diversity
3- Population genetic differentiation

Author: Anaisa Cajigas Gandia
License: MIT
Last updated: 8.4.2025

# R script

source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.10.txt")
getpackages() ## this is to get all required packages and dependencies sambaR need, and load them

setwd("path/to/directory")
list.files()

library(adegenet)
library(vcfR)
library (LEA)

# Import data
importdata(inputprefix="granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK",
           samplefile="populations.ingroup.tab", pop_order=c("DAM","ESQ","SAN","PAL"), 
           colourvector=c("green","springgreen4","coral","red4"),
           geofile="geofile.txt", do_citations=FALSE)

# Load sample sheet with populations information
samples<-read.csv("populations.ingroup.tab", sep="\t",header=T)

# Retain al data because it was previously filtered in GATK and PLINK)
filterdata(indmiss=1,snpmiss=1,min_mac=0,dohefilter=FALSE,min_spacing=0,
           nchroms=NULL,TsTvfilter=NULL)

#add2inds(samplefile="sample-info-geofile.txt",filteronly=FALSE, geomap_thres = -2)
inds$type<-ifelse(inds$pop=="DAM"|inds$pop=="ESQ",TRUE,FALSE)

findstructure(Kmin=1,Kmax=4,add_legend=TRUE, quickrun=FALSE,legend_pos="right",legend_cex=2,pop_order=NULL)
inferdemography(do_LEA=TRUE,Kmin=1,Kmax=4,do_f3=TRUE,f3_preparefiles=FALSE,jk_blocksize=1000)
LEAce(min_demes=1,max_demes=4,runanalysis=TRUE,nruns=50,export="pdf")
calcdiversity(nrsites=NULL, do_sfs=FALSE,do_venn=FALSE)
calcdistance(nchroms=NULL)
