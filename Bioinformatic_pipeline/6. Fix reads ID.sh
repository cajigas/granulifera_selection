### Fix reads ID ###
The script for both refrences is identical, only the input file and loaded modules change.

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 17.4.2025

#!/bin/bash
# Modules used for the genome
module load anaconda3/2019.03
module load java/16
conda init bash

# Modules used for the supertranscriptome
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load gatk/4.4.0.0

workDir="path/to/directory"
cd $workDir
for bamfile in $workDir/*_merged.sorted.deduped.bam ; do
    sample_name=$(basename -s _merged.sorted.deduped.bam $bamfile)
    echo -e "["$(date)"]\tRenaming.." $bamfile
gatk AddOrReplaceReadGroups -I $bamfile -O ${bamfile%%.bam}".ok.bam" \
       --RGID $sample_name \
       --RGLB PairedEnd \
       --RGPL Illumina \
       --RGPU merged \
       --RGSM $sample_name
done
