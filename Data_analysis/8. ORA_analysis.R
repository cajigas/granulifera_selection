### Over Representation Analysis ###

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 5.11.2025

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

install.packages("genekitr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library(genekitr)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(grid)
library (geneset)

# First step:Prepare input ID
genList <- read.csv("genome.csv", sep="\t",header=T)
gene_id <- genList$gene_name

# Remove sufixes to match standard gene symbols
clean_symbols <- gsub("(_[A-Z])$", "", gene_id)

# Transform gene list to symbols format on GO
symbol_id <- transId(clean_symbols, transTo =  "sym", keepNA = TRUE, unique = TRUE)  

# Second step: Get gene sets 

# bp: biological processes
# mf: molecular functions
# cc: celular components

gs_hum_bp <-  geneset::getGO(org = "human", ont = "bp")
gs_hum_mf <-  geneset::getGO(org = "human", ont = "mf")
gs_hum_cc <-  geneset::getGO(org = "human", ont = "cc")

# Third step: Perform ORA analysis
ora <- genORA(clean_symbols, geneset = gs_hum_cc)
ora1 <- genORA(clean_symbols, geneset = gs_hum_bp)
ora2 <- genORA(clean_symbols, geneset = gs_hum_mf)

# Fourth step: Plot ORA results_ genome

# Define a reusable theme
custom_theme <- theme(
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  axis.text.y = element_text(size = 14, margin = margin(b = 5)), # add margin for more spacing
  axis.text.x = element_text(size = 14),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title = element_text(size = 18, face = "bold")
)

# Individual enriched plots
p_cc <- plotEnrich(ora, plot_type = "bar", up_color = "red", down_color = "blue",
                   term_metric = "FoldEnrich", stats_metric = "pvalue") + 
  ggtitle("Cellular Component") + custom_theme

p_bp <- plotEnrich(ora1, plot_type = "bar", up_color = "red", down_color = "blue",
                   term_metric = "FoldEnrich", stats_metric = "pvalue") + 
  ggtitle("Biological Process") + custom_theme

p_mf <- plotEnrich(ora2, plot_type = "bar", up_color = "red", down_color = "blue",
                   term_metric = "FoldEnrich", stats_metric = "pvalue") + 
  ggtitle("Molecular Function") + custom_theme

# Combine vertically with more height
combined_genome <- p_mf / p_cc / p_bp + 
  plot_annotation(tag_levels = "A", title = "GO Enrichment across Ontologies_Genome") & 
  theme(plot.margin = margin(10, 10, 10, 10))

# Print with larger vertical space
print(combined_genome)
ggsave("combined_genome_plot.pdf", combined_genome, width = 10, height = 18) 

# Same steps were followed to conduct ORA analysis using the dataset obtained when
# mapping against the supertranscriptome 
