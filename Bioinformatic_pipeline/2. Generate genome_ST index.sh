### Generation of genome index ###

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

STAR --runThreadN 32 \
--runMode genomeGenerate --genomeDir \
--genomeFastaFiles path/to/reference_genome/GCA_033576555.1_ASM3357655v1_genomic.fa

### Generation of SuperTranscriptome index ###

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir \
--genomeFastaFiles path/to/reference_supertranscriptome/SuperDuper.sorted.fasta --genomeSAindexNbases 12 \
--sjdbGTFfile SuperDuperTrans.sorted.gff \
--sjdbOverhang 149 --limitGenomeGenerateRAM 100000000000
