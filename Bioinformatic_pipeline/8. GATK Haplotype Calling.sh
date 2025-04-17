### GATK Haplotype Calling ###
# The code for both references is the same, just modules change

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

# for the genome
#!/bin/bash
module load anaconda3/2019.03
module load java/8
module load gatk
conda init bash

samplesheet="path/to/file/granulifera_samples_merged_sorted_deduped.ok.cigar.bam.txt"
threads=$SLURM_JOB_CPUS_PER_NODE
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
BAM=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
reference="path/to/reference_genome/GCA_033576555.1_ASM3357655v1_genomic.fa" , or, "path/to/reference_supertranscriptome/SuperDuper.sorted.fasta"
INDIR="path/to/input_directory"
OUTDIR="path/to/output_directory"
cd $INDIR

gatk HaplotypeCaller --java-options '-Xmx160G' --native-pair-hmm-threads 1 \
-R $reference -I $BAM -ERC GVCF -O $OUTDIR/$name".vcf" --dont-use-soft-clipped-bases false 

# for the SuperTranscriptome
#!/bin/bash
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load gatk/4.4.0.0
