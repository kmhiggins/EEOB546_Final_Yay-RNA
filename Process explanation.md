To generate count tables used for downstream analysis, download the B73v5 genome fasta file, and gff3 annotation file from maizegdb.org. 
You will then want to convert the gff3 file to gtf using the following command from the Cufflinks package:

` gffread Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -T -o B73v5.gtf` 

This will convert the gff3 file into the gtf file required by STAR and featureCounts in the mRNA.sh and sRNA.sh scripts.

You will also want to generate the index files for Bowtie and STAR ahead of running the shell scripts.
These may take a long time to run.

Bowtie Index:

`bowtie-build Zm-B73-REFERENCE-NAM-5.0.fa Zea_mays`

STAR Index:

`STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /path/to/directory/B73_v5 --genomeFastaFiles /path/to/file/Zm-B73-REFERENCE-NAM-5.0.fa --sjdbGTFfile /path/to/file/B73v5.gtf --sjdbOverhang 99`

After generating the index files you'll want to run the shell scripts `count_table_mRNA.sh` and `count_table_sRNA.sh` on a high performance computing cluster using 

Once the count tables are generated through featureCounts, they are ready to be passed to R using Normalization_and_Heatmap.Rmd which uses DESeq2 to normalize and analyze the combined count tables. 

You will find the R script used for DESeq in the scripts folder as `Normalization_and_heatmap.Rmd`

the normalized files are saved under data as `normalized_mRNA.txt` and `normalized_sRNA.txt` and are required for the other R scripts included in this repository as are the files `mRNA_metadata.txt` and sRNA_metadata.txt` 

If you would like to run your own normalization on these, the raw counts are also included in the data folder. 
