### Merge tissues per individuals ###
# The code for both references is identical, just the modules change

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

# for the genome
#!/bin/bash
module load samtools

# for the SuperTranscriptome
#!/bin/bash
module load miniconda3/22.11.1  
module load gcc/13.2.0 
module load openmpi/4.1.6
module load curl/8.4.0-aixpq3w
module load samtools/1.17

cd path/to/input_directory
samtools merge --threads 8 DAM1_merged.sorted.deduped.bam DAM1*.sorted.deduped.bam
samtools merge --threads 8 DAM2_merged.sorted.deduped.bam DAM2*.sorted.deduped.bam
samtools merge --threads 8 DAM3_merged.sorted.deduped.bam DAM3*.sorted.deduped.bam
samtools merge --threads 8 DAM4_merged.sorted.deduped.bam DAM4*.sorted.deduped.bam
samtools merge --threads 8 DAM5_merged.sorted.deduped.bam DAM5*.sorted.deduped.bam
samtools merge --threads 8 DAM6_merged.sorted.deduped.bam DAM6*.sorted.deduped.bam
samtools merge --threads 8 DAM7_merged.sorted.deduped.bam DAM7*.sorted.deduped.bam
samtools merge --threads 8 DAM8_merged.sorted.deduped.bam DAM8*.sorted.deduped.bam
samtools merge --threads 8 ESQ1_merged.sorted.deduped.bam ESQ1*.sorted.deduped.bam
samtools merge --threads 8 ESQ2_merged.sorted.deduped.bam ESQ2*.sorted.deduped.bam
samtools merge --threads 8 ESQ3_merged.sorted.deduped.bam ESQ3*.sorted.deduped.bam
samtools merge --threads 8 ESQ4_merged.sorted.deduped.bam ESQ4*.sorted.deduped.bam
samtools merge --threads 8 ESQ5_merged.sorted.deduped.bam ESQ5*.sorted.deduped.bam
samtools merge --threads 8 ESQ6_merged.sorted.deduped.bam ESQ6*.sorted.deduped.bam
samtools merge --threads 8 ESQ7_merged.sorted.deduped.bam ESQ7*.sorted.deduped.bam
samtools merge --threads 8 ESQ8_merged.sorted.deduped.bam ESQ8*.sorted.deduped.bam
samtools merge --threads 8 PAL1_merged.sorted.deduped.bam PAL1*.sorted.deduped.bam
samtools merge --threads 8 PAL2_merged.sorted.deduped.bam PAL2*.sorted.deduped.bam
samtools merge --threads 8 PAL3_merged.sorted.deduped.bam PAL3*.sorted.deduped.bam
samtools merge --threads 8 PAL4_merged.sorted.deduped.bam PAL4*.sorted.deduped.bam
samtools merge --threads 8 PAL5_merged.sorted.deduped.bam PAL5*.sorted.deduped.bam
samtools merge --threads 8 PAL6_merged.sorted.deduped.bam PAL6*.sorted.deduped.bam
samtools merge --threads 8 PAL7_merged.sorted.deduped.bam PAL7*.sorted.deduped.bam
samtools merge --threads 8 PAL8_merged.sorted.deduped.bam PAL8*.sorted.deduped.bam
samtools merge --threads 8 SAN1_merged.sorted.deduped.bam SAN1*.sorted.deduped.bam
samtools merge --threads 8 SAN2_merged.sorted.deduped.bam SAN2*.sorted.deduped.bam
samtools merge --threads 8 SAN3_merged.sorted.deduped.bam SAN3*.sorted.deduped.bam
samtools merge --threads 8 SAN4_merged.sorted.deduped.bam SAN4*.sorted.deduped.bam
samtools merge --threads 8 SAN5_merged.sorted.deduped.bam SAN5*.sorted.deduped.bam
samtools merge --threads 8 SAN6_merged.sorted.deduped.bam SAN6*.sorted.deduped.bam
samtools merge --threads 8 SAN7_merged.sorted.deduped.bam SAN7*.sorted.deduped.bam
samtools merge --threads 8 SAN8_merged.sorted.deduped.bam SAN8*.sorted.deduped.bam

