


getwd()
setwd("C:/Users/Dasari Manasa/OneDrive/Documents/deseq2")

library(BiocManager)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(GEOquery)
library(tidyverse)

dir()
#step1 :preparing counts data

counts_data <- read.csv("raw_counts2.csv", row.names = 1)

# get sample info
geoid <- "GSE210029"

ges <- getGEO(geoid, GSEMatrix = TRUE)
phenodata <- pData(phenoData(ges[[1]]))

sample_info <- phenodata[ ,c(1,13)]

sample_info <- sample_info%>%
  rename(treatment = characteristics_ch1.3)%>%
  mutate(treatment =gsub( "treatment: ", "" , treatment))


# read the sample info



write.csv(sample_info,"sample_info.csv", row.names = FALSE)

sample_info <- read.csv("sample_info.csv",row.names = 1)

# making shore the row_name in sample_info matchs to coloum names in countsdata

all(colnames(counts_data) %in% rownames(sample_info))

# are they are in same order?
all(colnames(counts_data) == rownames(sample_info))


dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = sample_info,
                              design= ~ treatment)                                     
                                      
# prefiltaring : removing rows whith low gene count 

# keeping rows that have at least 10 read counts
keep <- rowSums(counts(dds))>=10

dds <- dds[keep,]
dim(dds)

# set factor level

dds$treatment <- relevel(dds$disease ,ref = "ntm_discase")
dds$treatment <- factor(dds$treatment)

# step3 run DESeq2

dds <- DESeq(dds)

deseq_results <- results(dds)

# changing deseq object R object data frame

deseq_results <- as.data.frame(deseq_results)
class(deseq_results)
head(deseq_results)

# order the results table by increasing p value

deseq_results_order <- deseq_results[order(deseq_results$pvalue),]

head(deseq_results_order)

# make some quaries

deseq_results["ENSG00000188171.10",]

# extract the most differetially expressed genes 
#select genes with a significant change in gene expression(adjusted p value below 0.05)
# and log 2 fold chang 1< and >1

# step1 filter based on p adjusted value

filtered <- deseq_results %>% filter(deseq_results$padj < 0.05) 

# filtered based on fold change , head we will use  threshold of 1

filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1 )

dim(deseq_results)
dim(filtered)

# save deseq results ,save the originaldata and filterd one

write.csv(deseq_results ,"deseq_results_all.csv")
write.csv(filtered , "deseq_results_filtered.csv")

# normalized counts

normalized_counts <- counts(dds ,normalize = TRUE)
head(normalized_counts)


write.csv(normalized_counts ,"normalized_counts.csv")

### visualization
# dispersion plot

plotDispEsts(dds)

# pca









                                      