### Mark duplicates ###
# This step was identical for both references, only the modules differ due to remote server software updates

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

# for genome
#!/bin/bash
module load anaconda3/2019.03
module load java/16
conda init bash

samplesheet="path/to/file/granulifera_samples_sorted.bam.txt"
threads=$SLURM_JOB_CPUS_PER_NODE
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
BAM=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
INDIR="path/to/input_directory"
OUTDIR="path/to/output_directory"
DEDUPED=$OUTDIR/"${sample}".sorted.deduped.bam
METRICS=$OUTDIR/"${sample}".picard-output.metrics

cd $OUTDIR
gatk MarkDuplicates -I $BAM -O $DEDUPED -VALIDATION_STRINGENCY SILENT -M $METRICS 

# for SuperTranscriptome
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load gatk/4.4.0.0

threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
gatk MarkDuplicates -I $BAM -O $DEDUPED -VALIDATION_STRINGENCY SILENT -M $METRICS 



