### Generation of genome index ###

Author: Anaisa Cajigas Gandia
License: MIT
Last updated: 8.4.2025

STAR --runThreadN 32 \
--runMode genomeGenerate --genomeDir \
--genomeFastaFiles GCA_033576555.1_ASM3357655v1_genomic.fa


### Generation of SuperTranscriptome index ###

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir \
--genomeFastaFiles SuperDuper.sorted.fasta --genomeSAindexNbases 12 \
--sjdbGTFfile SuperDuperTrans.sorted.gff \
--sjdbOverhang 149 --limitGenomeGenerateRAM 100000000000
