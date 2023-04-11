#!/bin/bash
SECONDS=0

## step 1 ###
fastqc raw_data/SRR20853616/*.fastq -o raw_data/

#runtrimmamatic ##

trimmomatic PE -threads 4 raw_data/SRR20853616/SRR20853616_1.fastq  raw_data/SRR20853616/SRR20853616_2.fastq trimmdata/SRR20853616_1_paired.fastq trimmdata/SRR20853616_1_unpaired.fastq trimmdata/SRR20853616_2_paired.fastq trimmdata/SRR20853616_2_unpaired.fastq ILLUMINACLIP:/trimmdata/TruSeq3-PE-2.fa:2:30:10 HEADCROP:10 MINLEN:36 SLIDINGWINDOW:4:15

#fastqc trimmdata/SRR20853616_1_paired.fastq 
#fastqc trimmdata/SRR20853616_2_paired.fastq

#step2 ##

hisat2 -x hisat2/grch38/genome -1 trimmdata/SRR20853616_1_paired.fastq -2 trimmdata/SRR20853616_2_paired.fastq | samtools sort -o hisat2/demo.trimmed.bam

# Step - 3 run -feature count quantification ###
# get gtf

#wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

featureCounts -T 4 -p -a hisat2/Homo_sapiens.GRCh38.109.gtf -o quants/featurecounts.txt hisat2/demo.trimmed.bam

duration=$SECONDS
echo "$(($duration/60)) minutes and $(($duration % 60 )) SECONDS. "
