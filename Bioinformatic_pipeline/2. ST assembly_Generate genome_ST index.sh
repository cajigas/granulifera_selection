### Generation of ST assembly, and ST and genome indexes ###

Author: Anaisa Cajigas Gandia
License: GNU
Last updated: 8.4.2025

### Generation of genome index ###
STAR --runThreadN 32 \
--runMode genomeGenerate --genomeDir \
--genomeFastaFiles path/to/reference_genome/GCA_033576555.1_ASM3357655v1_genomic.fa

### SuperTranscriptome assembly with Lace ###

/data/Software/necklace-1.01/tools/bin/lace transcript_collection.okay.mrna.annotations-no-contams.shortheaders.fasta
/data/granulifera/all_tissues/LACE/annotations-no-contams_trancript-to-gene_mappings.tsv \
--outputDir /data/granulifera/all_tissues/LACE/output
--cores 24 --alternate --maxTran 500

### Generation of SuperTranscriptome index ###

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir \
--genomeFastaFiles path/to/reference_supertranscriptome/SuperDuper.sorted.fasta --genomeSAindexNbases 12 \
--sjdbGTFfile SuperDuperTrans.sorted.gff \
--sjdbOverhang 149 --limitGenomeGenerateRAM 100000000000
