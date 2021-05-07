# EEOB546_Final_Yay-RNA

In this repository you will find shell scripts, R scripts, and data involved in the attempted duplication of the paper "Variation and Inheritance of Small RNAs in Maize
Inbreds and F1 Hybrids" by Crisp, et al. 2020. 

## Required Packages
The following shell script packages are required to generate count tables:

`sra-toolkit`

`trimmomatic`

`STAR`

`subread`

`bowtie`

The following R packages are required to generate count tables and heatmaps:
> purrr
> 
> tidyverse
> 
> GenomicRanges
> 
> Biobase
> 
> DESeq2
> 
> limma
> 
> Glimma
> 
> gplots
> 
> RColorBrewer
> 
> GEOquery
> 
> tibble
> 
> ggplot2
> 
> reshape2
> 
> viridis

## Data Folder

In this folder you will find data required for several parts of processing

#### Getting counts

You will need the following data files for generating read counts:

> mRNA_accessions.txt
> 
> sRNA_Accessions.txt

#### R Studio
You will need the following data files for running Rscripts:

> mRNA_metadata.txt
> 
> sRNA_metadata.txt

If you choose to generate your own count files the following are not necessary, but are if you are only looking to run the R scripts:

> mRNA_raw_counts.txt
> 
> sRNA_raw_counts.txt
> 
> normalized_mRNA.txt
> 
> normalized_sRNA.txt

## Scripts Folder

In this folder you will find two types of scripts used in our project, shell scripts and R scripts

The shell scripts, count_table_mRNA.sh and count_table_sRNA.sh, are used for the initial mapping step and would be best run on a an HPC. They are dependent on the availability of data files mRNA_accessions.txt and sRNA_accessions.txt

The R script normalization_and_heatmap.Rmd is used for normalizing the counts, running them through DESeq2 and generating heatmaps of differentially expressed genes. It is dependent on either the featureCount.txt files generated from the shell scripts or on the normalized or raw_count data. 



The R scripts PCA_of_mRNA.Rmd and ANOVA_of_sRNA.Rmd were used to attempt to recreate figures 3 and 4 in Crisp et al., 2020. 

# Required Packages
The following packages are needed to create the PCA:

>BiocManager
>
>pcaMethods 

You will also need the following files found in data:

>normalized_mRNA.txt 

The following packages are needed to create the ANOVA: 

>ggpubr
>
>car
>
>dplyr
>
>tidyverse 
>
>stringr 

You will also need the following files found in data:

>normalized_mRNA.txt
>
>mRNA_metadate.txt




