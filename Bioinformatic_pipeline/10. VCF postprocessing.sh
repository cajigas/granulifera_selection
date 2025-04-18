### VCF postprocessing with vcftools and PLINK ###
# The code is identical for both references, just file names change

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

help> https://www.cog-genomics.org/plink/2.0/data#make_pgen
help> https://github.com/vcftools/vcftools 

# Filter in vcftools

vcftools --vcf granulifera_GATK_variants.filtered.SNPS.biallelic.vcf --minDP 3 --maf 0.05 --max-missing 0.75 \
--recode --recode-INFO-all --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8

# Remove loci in Linkage disequilibrium (LD) with PLINK
plink2 --vcf granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.recode.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --indep-pairwise 50 1 0.7 --bad-ld 

plink2 --vcf granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.recode.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --exclude plink2.prune.out --make-bed \
--recode --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK

# Export in different formats
plink --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK \
--allow-extra-chr --recode A --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK

plink --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8_noCDS.LDpruned.PLINK \
--allow-extra-chr --recode --tab --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK

plink2 --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK --allow-extra-chr --freq
plink2 --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK --allow-extra-chr --make-rel square --read-freq plink2.afreq

plink2 --bfile granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK \
--allow-extra-chr --export vcf --out granulifera_GATK_variants.filtered.biallelic.SNPs.minDP3mac3maxmissing8.LDpruned.PLINK 
