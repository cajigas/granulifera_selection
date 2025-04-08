### Fix groups ID ###
# The code for both references is identical, just modules change

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

# for the genome
#!/bin/bash
module load anaconda3/2019.03
module load java/16
conda init bash

# for the SuperTranscriptome
#!/bin/bash
module load miniconda3/22.11.1
module load gcc/13.2.0
module load openmpi/4.1.6
module load openjdk/17.0.8.1_1
module load gatk/4.4.0.0


cd $workDir
for bamfile in workdir/*_merged.sorted.deduped.bam ; do
    sample_name=$(basename -s _merged.sorted.deduped.bam $bamfile)
    echo -e "["$(date)"]\tRenaming.." $bamfile
gatk AddOrReplaceReadGroups -I $bamfile -O ${bamfile%%.bam}".ok.bam" \
       --RGID $sample_name \
       --RGLB PairedEnd \
       --RGPL Illumina \
       --RGPU merged \
       --RGSM $sample_name
done
