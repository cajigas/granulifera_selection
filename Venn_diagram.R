### Venn diagram ###
This script shows how the Venn diagrams to represent common SNPs across the three selection scans for the genome and the superTranscriptome.
The script for both references is identical, just the input file changes.



setwd("path/to /directory")
library(qqplot)
library(VennDiagram)

# Load file with SNPs identified in the genome
venn <- read.csv("venn_diagram.csv",header=TRUE)
# Load file with SNPs identified in the superTranscriptome
venn <- read.csv("Venn_diagram_ST.csv" ,header=TRUE)

# Identify unique SNPs per selection scan
set1 <- unique(venn$SNP_Ohana)
set2 <- unique(venn$SNP_Bayescan)
set3 <- unique(venn$SNP_PCAdapt)

# Remove empty strings and NA values
set1 <- set1[trimws(set1) != "" & !is.na(set1)]
set2 <- set2[trimws(set2) != "" & !is.na(set2)]
set3 <- set3[trimws(set3) != "" & !is.na(set3)]

# Plot the Venn diagram
venn.plot <- venn.diagram(     # or Venn_diagram_ST
  x = list(A = set1, B = set2, C = set3),
  category.names = c("Ohana", "BayeScan", "PCAdapt"),
  filename = NULL,      
  fill = c("blue", "red", "green"), 
  alpha = 0.5,           
  cex = 2.5,             
  cat.cex = 2,         
  cat.pos = c(-20, 20, 180) 
)
grid.draw(venn.plot)

