#!/bin/bash
#sbatch --array=0-104 count_table_mRNA.sh
#SBATCH --time=72:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2 
#SBATCH --mem=80G
#SBATCH --job-name="Yay-RNA" 
#SBATCH --mail-user=kaitlinj@iastate.edu   # email address
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL

#cd /work/LAS/sna-lab/kmh/BCB585_tmp
#Load modules
module load sra-toolkit/2.9.6-ub7kz5h
module load py-cutadapt/1.13-py2-e2j6srt #apparently necessary to run readarray command
module load trimmomatic/0.36-lkktrba
module load star/2.5.3a-rdzqclb
module load r-rjava/0.9-8-py2-r3.4-wadatwr
module load samtools
module load subread/1.6.0-ak6vxhs

# Read in sample file to use parts
readarray -t FILES < mRNA_accessions.txt
SpudsMac=${FILES[$SLURM_ARRAY_TASK_ID]}
IFS=' '
read -ra INFO <<< "$SpudsMac"
# Cut relevant column information into variables
SRA="$(cut -f 1 <<< $INFO)"
SAMPLE="$(cut -f 2 <<< $INFO)"

# Pull FASTA information from NCBI based on SRA number, split-files will split the files into 3 with two of them being _1 and _2 which is useful for trimmomatic
fasterq-dump -p -t /tmp --split-files ${SRA}
# Trim adapters using trimmomatic
trimmomatic PE ${SRA}_1.fastq ${SRA}_2.fastq ${SAMPLE}_out_1P.fq ${SAMPLE}_out_1U.fq ${SAMPLE}_out_2P.fq ${SAMPLE}_out_2U.fq ILLUMINACLIP:TruSeq3-PE.fa:4:30:10 MINLEN:30

# Run mapping
STAR --runMode alignReads --runThreadN 6 --genomeDir /B73_v5 --readFilesIn ${SAMPLE}_out_1P.fq ${SAMPLE}_out_2P.fq --twopassMode Basic --outFileNamePrefix ${SAMPLE}

# Run HTSeq to generate unique count tables
featureCounts -p -t exon -g gene_id -a Zea_mays.AGPv4.33.gtf -o ${SAMPLE}_featureCounts.txt ${SAMPLE}Aligned.out.sam

