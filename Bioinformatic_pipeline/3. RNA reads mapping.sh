### Mapping reads against reference genome using STAR ###

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

#!/bin/bash
module load gcc/9.3.0
module load openmpi/gcc.9/3.1.5

threads=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS

samplesheet="path/to/file/granulifera_samples_trimmed.txt"
samplename=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`
genomeDir="path/to/reference_genome"
workDir= "path/to/working_directory"
mkdir $workDir
cd $workDir

STAR --readFilesCommand zcat --outFileNamePrefix $samplename \
	--outSAMtype BAM Unsorted --outSAMmapqUnique 60 \
	--outSAMattrRGline ID:$samplename CN:ZoolInst_TiHo \
	LB:PairedEnd PL:Illumina PU:Unknown SM:$samplename \
	--genomeDir $genomeDir --runThreadN $SLURM_NTASKS --readFilesIn $r1 $r2 --twopassMode Basic 
	
samtools sort -@ 2 -o $samplename"_sorted.bam" $samplename"Aligned.out.bam"
samtools index -c -@ 2 $samplename"_sorted.bam"
rm $workDir/$samplename"Aligned.out.bam"

### Mapping reads against reference SuperTranscriptome using STAR ###

# Note that by the time when the analysis with the SuperTranscriptome was done, there were software updates in the supercomputer remote server
# and therefore, new modules and module versions had to be used

#!/bin/bash
module load miniconda3/22.11.1
module load openjdk/17.0.8.1_1
module load openmpi/4.1.6
module load gcc/13.2.0
module load curl/8.4.0-aixpq3w
module load samtools/1.17

threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

samplesheet="path/to/file/granulifera_samples_trimmed.txt"
samplename=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`
supertranscriptDir="path/to/reference_supertranscriptome"
workDir= "path/to/working_directory"
cd $workDir

STAR --readFilesCommand zcat --outFileNamePrefix $samplename \
--outSAMtype BAM Unsorted --outSAMmapqUnique 60 --outSAMstrandField intronMotif \
--outSAMattrRGline ID:$samplename CN:ZoolInst_TiHo LB:PairedEnd PL:Illumina PU:Unknown SM:$samplename \
--genomeDir $supertranscriptDir --runThreadN $SLURM_CPUS_PER_TASK --readFilesIn $r1 $r2 --twopassMode Basic \
--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5	

samtools sort -@ $SLURM_CPUS_PER_TASK -o $samplename"_sorted.bam" $samplename"Aligned.out.bam"
samtools index -c -@ $SLURM_CPUS_PER_TASK $samplename"_sorted.bam"
rm $workDir/$samplename"Aligned.out.bam"

