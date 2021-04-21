#!/bin/bash
#sbatch --array=0-56 count_table_sRNA.sh
#SBATCH --time=72:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
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
module load gcc/7.3.0-xegsmw4
module load bowtie/1.2-26wthei
module load r-rjava/0.9-8-py2-r3.4-wadatwr
module load samtools
module load py-htseq/0.11.2-py2-4757mqt
module load subread/1.6.0-ak6vxhs
# Read in sample file to use parts
readarray -t FILES < sRNA_Accessions.txt
SpudsMac=${FILES[$SLURM_ARRAY_TASK_ID]}
IFS=' '
read -ra INFO <<< "$SpudsMac"
# Cut relevant column information into variables
SRA="$(cut -f 1 <<< $INFO)"
SAMPLE="$(cut -f 2 <<< $INFO)"
# Pull FASTA information from NCBI based on SRA number
fasterq-dump -p -t /tmp ${SRA}
# Trim adapters using trimmomatic
trimmomatic SE -phred33 ${SRA}.fastq ${SAMPLE}_1_out.fastq ILLUMINACLIP:TruSeq3-SE.fa:4:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
#Align sRNAs using bowtie
bowtie -v 1 -m 50 Zea_mays ${SAMPLE}_1_out.fastq -S ${SAMPLE}.sam
#Generate count tables
featureCounts -M -d 18 -D 34 -t exon -g gene_id -a B73v5.gtf -o ${SAMPLE}_featureCounts.txt ${SAMPLE}.sam