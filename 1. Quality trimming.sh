### RNA reads quality trimming ###

Author: Anaisa Cajigas Gandia
License: MIT
Last updated: 8.4.2025

#!/bin/bash

module load openmpi/gcc.9/3.1.5

# provide the directory with the reads and 
# a tab separated table of PE reads samples in three columns (id mate-1 mate-2)
READS_DIR= "my directory"
samplesheet="granulifera_samples_raw.txt"
threads=$SLURM_NTASKS
samplename=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`

#fastp array (trimming many sequences in a same directory)
fastp --thread $SLURM_NTASKS --detect_adapter_for_pe \
--in1 $r1 --in2 $r2 \
--out1 ${r1%%.fastq.gz}".trim.fastq.gz" \
--out2 ${r2%%.fastq.gz}".trim.fastq.gz" \
--html ${samplename}".fastp.html" \
--json ${samplename}".fastp.json" 

