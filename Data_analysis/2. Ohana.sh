### Selection scan wih Ohana ###
The script for both references is identical, just the input file changes.
Ohana scans are done in a Linux server, and the results are analysed in R.

Author: Ariel Rodriguez
License: GNU
Last updated: 17.4.2025

# Running Ohana in Linux server #
help: https://github.com/jade-cheng/ohana/wiki/Population-or-ancestry-specific-selection-scan
 
# Installing Ohana
cd Software
sudo apt install libblas-dev liblapacke-dev
git clone https://github.com/jade-cheng/ohana
cd ./ohana
make

mkdir "path/to/output_directory"
cd "path/to/output_directory"

plink --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8_noCDS.LDpruned.PLINK \
--allow-extra-chr --recode12 --geno 0.0 --tab --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing_noCDS.LDpruned.PLINK

ohana/bin/convert ped2dgm granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing_noCDS.LDpruned.PLINK.ped granulifera.LDpruned.dgm
ohana/bin/qpas granulifera.LDpruned.dgm -k 3 -qo q_pruned.matrix -fo f_pruned.matrix -mi 500
ohana/bin/nemeco granulifera.LDpruned.dgm f_pruned.matrix -mi 500 -co count_pruned.matrix
ohana/bin/convert cov2nwk count_pruned.matrix granulifera_k3.tre

# Plot
python Software/ohana/tools/plot-q.py q_pruned.matrix q-bar-chart_k3.pdf

# Now we do the qpas for the full dataset using as input the results from unlinked loci
plink --vcf granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.recode.vcf \
--allow-extra-chr --recode12 --geno 0.0 --tab --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing.OHANA

ohana/bin/convert ped2dgm granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing.OHANA.ped granulifera.FULL.dgm
ohana/bin/qpas granulifera.FULL.dgm -k 3 -qi q_pruned.matrix -fo f_full.matrix -mi 500 -fq -e 0.0001
ohana/bin/selscan granulifera.FULL.dgm f_full.matrix count_pruned.matrix > granulifera_ohana_selscanK3_all.tsv

# Now we do population-specific scans
ohana/bin/selscan granulifera.FULL.dgm f_full.matrix count_pruned.matrix -cs cs_matrix_k0 > granulifera_ohana_selscanK3_greens.tsv
paste granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing.OHANA.snpids granulifera_ohana_selscanK3_greens.tsv > granulifera_ohana_selscanK3_greens.snps.tsv

ohana/bin/selscan granulifera.FULL.dgm f_full.matrix count_pruned.matrix -cs cs_matrix_k1 > granulifera_ohana_selscanK3_PAL.tsv
paste granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing.OHANA.snpids granulifera_ohana_selscanK3_PAL.tsv > granulifera_ohana_selscanK3_PAL.snps.tsv

ohana/bin/selscan granulifera.FULL.dgm f_full.matrix count_pruned.matrix -cs cs_matrix_k2 > granulifera_ohana_selscanK3_SAN.tsv
paste granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3nomissing.OHANA.snpids granulifera_ohana_selscanK3_SAN.tsv > granulifera_ohana_selscanK3_SAN.snps.tsv

# After this step, for the superTranscriptome data, we do:
# combine with SNPids

plink2 --vcf granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.recode.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --make-bed --recode --geno 0.0 --out granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3NoMissing.SNPids

awk '{print $2}' granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3NoMissing.SNPids.bim \
> granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3NoMissing.SNPids.txt
paste granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.minDP3mac3NoMissing.SNPids.txt granuliferaST_ohana_selscanK3_greens.tsv \
> granuliferaST_ohana_selscanK3_greens.snps.tsv

# Ohana results post-processing in R #
setwd("path/to/directory")

library(qqman)
library(qvalue)

# Load data for green populations
ohana.greens<-read.csv("granulifera_ohana_selscanK3_greens.snps.tsv", sep="\t",row.names=NULL)
ohana.greens$llk<-log(ohana.greens$lle.ratio)
quant.green<-quantile(ohana.greens$llk,probs=0.99,na.rm=TRUE)
summary(ohana.greens[ohana.greens$llk>quant.green,])
write.csv(ohana.greens[ohana.greens$llk>quant.green,],
"granulifera_ohana_selscanK3_greens.TOP1.csv" )

# Load data for PAL population
ohana.PALs<-read.csv("granulifera_ohana_selscanK3_PAL.snps.tsv", sep="\t",row.names=NULL)
ohana.PALs$llk<-log(ohana.PALs$lle.ratio)
quant.PAL<-quantile(ohana.PALs$llk,probs=0.99,na.rm=TRUE)
summary(ohana.PALs[ohana.PALs$llk>quant.PAL,])
write.csv(ohana.PALs[ohana.PALs$llk>quant.PAL,],
"granulifera_ohana_selscanK3_PALs.TOP1.csv" )

# Load data for SAN population
ohana.SANs<-read.csv("granulifera_ohana_selscanK3_SAN.snps.tsv", sep="\t",row.names=NULL)
ohana.SANs$llk<-log(ohana.SANs$lle.ratio)
quant.SAN<-quantile(ohana.SANs$llk,probs=0.99,na.rm=TRUE)
summary(ohana.SANs[ohana.SANs$llk>quant.SAN,])
write.csv(ohana.SANs[ohana.SANs$llk>quant.SAN,],
"granulifera_ohana_selscanK3_SANs.TOP1.csv" )



