### SplitNCigar ###
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

samplesheet="path/to/file/granulifera_samples_merged_sorted_deduped.ok.bam.txt"
threads=$SLURM_JOB_CPUS_PER_NODE
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
BAM=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
reference="path/to/reference_genome/GCA_033576555.1_ASM3357655v1_genomic.fa"
INDIR="path/to/input_directory"
OUTDIR="path/to/output_directory"
DEDUPED=$OUTDIR/"${name}"_merged.sorted.deduped.ok.bam
CIGAR=$OUTDIR/"${name}"_merged.sorted.deduped.ok.cigar.bam

# Create sequence dictionary
samtools faidx -c $reference 
java -jar  picard.jar CreateSequenceDictionary R=GCA_033576555.1_ASM3357655v1_genomic.fa \
O= GCA_033576555.1_ASM3357655v1_genomic.dict

#Create directory for temporal big files that are generated during this step
mkdir temp
cd $OUTDIR
gatk SplitNCigarReads --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' \
--input $DEDUPED --output $CIGAR --reference $reference --skip-mapping-quality-transform \
--tmp-dir temp --create-output-bam-index false

# Create csi index
samtools index -c -@ 4 $CIGAR

# for the SuperTranscriptome
#!/bin/bash
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load curl/8.4.0-aixpq3w
module load samtools/1.17
module load gatk/4.4.0.0

# Create sequence dictionary
samtools faidx -c SuperDuper.sorted.fasta
java -jar picard.jar CreateSequenceDictionary R= SuperDuper.sorted.fasta \
O= SuperDuper.sorted.dict




