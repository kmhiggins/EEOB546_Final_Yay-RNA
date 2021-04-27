# EEOB546_Final_Yay-RNA

In this repository you will find shell scripts, R scripts, and some data involved in the attempted duplication of the paper "Variation and Inheritance of Small RNAs in Maize
Inbreds and F1 Hybrids" by Crisp, et al. 2020. 

## Required Packages
The following shell script packages are required to generate count tables:

`sra-toolkit`

`cutadapt`

`trimmomatic`

`STAR`

`subread`

`bowtie`

## Data 
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

Once the count tables are generated, they are ready to be passed to R using {R Script here}
