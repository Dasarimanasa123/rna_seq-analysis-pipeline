

# !/bash/bin
cat raw_data/SRR_Acc_List.txt | parallel "fastq-dump --gzip --split-files {}"

# fast qc

fastqc -t 64 raw_data/*.fastq.gz -o raw_data/

### multiqc ###

multiqc .

## step 2 - trimmomatic ###

cat raw_data/SRR_Acc_List.txt | parallel -j 4 "trimmomatic PE -threads 4 raw_data/{}_1.fastq.gz  raw_data/{}_2.fastq.gz  trimm_data/{}_1_paired.fastq.gz trimm_data/{}_1_unpaired.fastq.gz trimm_data/{}_2_paired.fastq.gz trimm_data/{}_2_unpaired.fastq.gz HEADCROP:10 MINLEN:36 SLIDINGWINDOW:4:15"


#step 3 ##

cat raw_data/SRR_Acc_List.txt | parallel -j 4 "hisat2 -x hisat2/grch38/genome -1 trimm_data/{}_1_paired.fastq.gz -2 trimm_data/{}_2_paired.fastq.gz | samtools sort -o hisat2/demo.trimmed.bam"

# Step - 3 run -feature count quantification ###

featureCounts -T 4 -p -a hisat2/Homo_sapiens.GRCh38.109.gtf -o quants/featurecounts.txt hisat2/demo.trimmed.bam

echo "rna seq completed"
