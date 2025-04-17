### GATK Merge individuals ###
# The code is the same for both references, just modules, reference and file name change accordingly

Author: Anaisa Cajigas Gandia
License: MIT
Last updated: 8.4.2025

# for the genome
gatk CombineGVCFs --java-options '-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=28' \
-R GCA_033576555.1_ASM3357655v1_genomic.fa \
--variant DAM1.vcf \
--variant DAM2.vcf \
--variant DAM3.vcf \
--variant DAM4.vcf \
--variant DAM5.vcf \
--variant DAM6.vcf \
--variant DAM7.vcf \
--variant DAM8.vcf \
--variant ESQ1.vcf \
--variant ESQ2.vcf \
--variant ESQ3.vcf \
--variant ESQ4.vcf \
--variant ESQ5.vcf \
--variant ESQ6.vcf \
--variant ESQ7.vcf \
--variant ESQ8.vcf \
--variant PAL1.vcf \
--variant PAL2.vcf \
--variant PAL3.vcf \
--variant PAL4.vcf \
--variant PAL5.vcf \
--variant PAL6.vcf \
--variant PAL7.vcf \
--variant PAL8.vcf \
--variant SAN1.vcf \
--variant SAN2.vcf \
--variant SAN3.vcf \
--variant SAN4.vcf \
--variant SAN5.vcf \
--variant SAN6.vcf \
--variant SAN7.vcf \
--variant SAN8.vcf \
-O granulifera_GATK_variants.gvcf

# Convert GVCF to VCF
gatk GenotypeGVCFs --java-options '-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=28' \
-R GCA_033576555.1_ASM3357655v1_genomic.fa \
-V granulifera_GATK_variants.gvcf -O granulifera_GATK_variants.vcf

# Variant Filtration
java -jar gatk-package-4.1.2.0-local.jar VariantFiltration -R GCA_033576555.1_ASM3357655v1_genomic.fa \
-V granulifera_GATK_variants.vcf -window 35 -cluster 3 --filter-name FSQD -filter "FS > 30.0 || QD < 2.0" \
-O  granulifera_GATK_variants.filtered.vcf 2> granulifera_GATK_variants.VariantFilter.log

java -Xmx60g -Xms60g -XX:+UseParallelGC -XX:ParallelGCThreads=28 -jar gatk-package-4.1.2.0-local.jar SelectVariants \
-R GCA_033576555.1_ASM3357655v1_genomic.fa -V granulifera_GATK_variants.filtered.vcf \
--restrict-alleles-to BIALLELIC --exclude-filtered true -O granulifera_GATK_variants.filtered.biallelic.vcf 2>granulifera_GATK_variants.filtered.SelectVariants_BIALLELIC.log

# Divide into SNPs and Indels
java -Xmx60g -Xms60g -XX:+UseParallelGC -jar gatk-package-4.1.2.0-local.jar SelectVariants \
-R GCA_033576555.1_ASM3357655v1_genomic.fa -V granulifera_GATK_variants.filtered.biallelic.vcf \
--select-type-to-include SNP -O granulifera_GATK_variants.filtered.biallelic.SNPs.vcf 

java -Xmx60g -Xms60g -XX:+UseParallelGC -jar gatk-package-4.1.2.0-local.jar SelectVariants \
-R GCA_033576555.1_ASM3357655v1_genomic.fa -V granulifera_GATK_variants.filtered.biallelic.vcf \
--select-type-to-include INDEL -O granulifera_GATK_variants.filtered.biallelic.INDELS.vcf
 
# for the SuperTranscriptome

module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load python/3.11.6
module load gatk/4.4.0.0

gatk CombineGVCFs --java-options '-Xmx700G -XX:+UseParallelGC -XX:ParallelGCThreads=96' \
-R SuperDuper.sorted.fasta \
--variant DAM1.vcf \
--variant DAM2.vcf \
--variant DAM3.vcf \
--variant DAM4.vcf \
--variant DAM5.vcf \
--variant DAM6.vcf \
--variant DAM7.vcf \
--variant DAM8.vcf \
--variant ESQ1.vcf \
--variant ESQ2.vcf \
--variant ESQ3.vcf \
--variant ESQ4.vcf \
--variant ESQ5.vcf \
--variant ESQ6.vcf \
--variant ESQ7.vcf \
--variant ESQ8.vcf \
--variant PAL1.vcf \
--variant PAL2.vcf \
--variant PAL3.vcf \
--variant PAL4.vcf \
--variant PAL5.vcf \
--variant PAL6.vcf \
--variant PAL7.vcf \
--variant PAL8.vcf \
--variant SAN1.vcf \
--variant SAN2.vcf \
--variant SAN3.vcf \
--variant SAN4.vcf \
--variant SAN5.vcf \
--variant SAN6.vcf \
--variant SAN7.vcf \
--variant SAN8.vcf \
-O granulifera_Supertranscriptome_variants.gvcf

#Convert GVCF to VCF
gatk GenotypeGVCFs --java-options '-Xmx700g -XX:+UseParallelGC -XX:ParallelGCThreads=96' \
-R SuperDuper.sorted.fasta \
-V granulifera_Supertranscriptome_variants.gvcf \
-O granulifera_Supertranscriptome_variants.vcf

#Variant filtration
gatk VariantFiltration -R SuperDuper.sorted.fasta \
-V granulifera_Supertranscriptome_variants.vcf --filter-name FSQD -filter "FS > 30.0 || QD < 2.0" \
-O granulifera_Supertranscriptome_variants.filtered.vcf 2> granulifera_Supertranscriptome_variants.VariantFilter.log

gatk SelectVariants --java-options '-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=32' \
-R SuperDuper.sorted.fasta \
-V granulifera_Supertranscriptome_variants.filtered.vcf --restrict-alleles-to BIALLELIC --exclude-filtered true \
-O granulifera_Supertranscriptome_variants.filtered.biallelic.vcf 2> granulifera_Supertranscriptome_variants_filtered.SelectVariants_BIALLELIC.log

#Divide into SNPs and Indels
gatk SelectVariants --java-options '-Xmx100g -Xms100g -XX:+UseParallelGC -XX:ParallelGCThreads=32' \
-R SuperDuper.sorted.fasta \
-V granulifera_Supertranscriptome_variants.filtered.biallelic.vcf --select-type-to-include SNP \
-O granulifera_Supertranscriptome_variants.filtered.biallelic.SNPs.vcf 

gatk SelectVariants --java-options '-Xmx100g -Xms100g -XX:+UseParallelGC -XX:ParallelGCThreads=32' \
-R SuperDuper.sorted.fasta \
-V granulifera_Supertranscriptome_variants.filtered.biallelic.vcf --select-type-to-include INDEL \
-O granulifera_Supertranscriptome_variants.filtered.biallelic.INDELS.vcf
