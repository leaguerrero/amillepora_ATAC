## First round data exploration

# Download packages:
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)

# File Set-up
file.dir <- "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Gene_count_files/"
my.files <- grep(".txt", list.files(file.dir), value=TRUE)
my.metadata <- read.table("/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)

# Create sample table
sampleNames <- c("01_27_Y","01_28_R","01_40_Y","01_42_R","01_45_Y","01_46_Y",
                 "01_51_Y","01_52_Y","01_54_Y","01_58_Y","01_63_Y","01_68_Y",
                 "01_78_R","01_92_Y","02_29_Y","02_34_Y","02_35_R","02_47_Y",
                 "02_65_R","02_79_R","03_29_R","03_48_Y","03_62_R","03_72_R",
                 "03_90_Y","03_91_R","04_48_R","04_65_Y","04_70_Y","04_95_Y",
                 "05_43_R","05_46_R","05_50_Y","05_62_Y","05_82_Y","05_88_R",
                 "06_22_Y","06_32_R","06_59_R","06_64_R","06_67_R","06_72_Y",
                 "06_82_R","06_92_R","07_31_R","07_36_R","07_39_R","07_53_Y",
                 "07_64_Y","07_88_Y","08_19_R","08_30_R","08_33_Y","08_36_Y",
                 "08_45_R","08_50_R","08_85_Y","08_95_R")

my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files, 
                             condition = my.metadata[,2:9])

# Convert variables to factors
factorVars <- c("condition.acclimation.txt", "condition.acclimation.group", 
                "condition.stress.txt", "condition.stress.group", 
                "condition.stress.day", "condition.ext.date",
                "condition.ext.batch", "condition.genotype")

my.sampleTable[, factorVars] <- lapply(my.sampleTable[, factorVars], factor)
colnames(my.sampleTable)

# Check my.sampleTable
my.sampleTable


# Create DESeqDataSet Object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.txt +
    condition.stress.txt +
    condition.stress.day +
    condition.ext.batch +
    condition.acclimation.txt:condition.stress.txt)

# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Differential Expression Analysis
dds <- DESeq(ddsHTSeq) 
resultsNames(dds)
# saveRDS(dds, file = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Data_out/dds.RData")

# Report: Number of transcripts
num_transcripts <- nrow(dds)
num_transcripts # 19649 on September 08 2023

# Save for topGO analysis
total.set.gene.names <- rownames(assay(dds))
total.set.gene.names.df <- data.frame(gene = total.set.gene.names)
# write.csv(total.set.gene.names.df, file = '~/Lab Notebook/Chapter2/Data_analysis/TagSeq/Data_out/total.set.gene.names.csv', row.names = FALSE)

# Variance stabilizing transformation
vtd <- vst(dds)

# PCA: heat stress day
deseq_PCA_heat_stress_day <- plotPCA(vtd,intgroup=c( "condition.stress.day"), returnData = TRUE) 
deseq_PCA_heat_stress_day$condition.stress.day

# Plot PCA Heat Stress day
PCA_plot_hs_day <- ggplot(deseq_PCA_heat_stress_day,aes(x=PC1,y=PC2, color=group)) +
  theme_classic(base_size = 20) +
  labs(title="Gene Expression") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs_day

# PCA: extraction batch
deseq_PCA_batch <- plotPCA(vtd,intgroup=c( "condition.ext.batch"), returnData = TRUE) 
deseq_PCA_batch$group<- as.factor(deseq_PCA_batch$group)

# Plot PCA ext batch
PCA_plot_ext_batch <- ggplot(deseq_PCA_batch,aes(x=PC1,y=PC2, color = group)) +
  theme_classic(base_size = 20) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="Gene Expression", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_ext_batch

# PCA: Heat stress treatment
deseq_PCA_hs <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs$group<- as.factor(deseq_PCA_hs$group)

# Plot PCA heat stress treatment
PCA_plot_hs <- ggplot(deseq_PCA_hs,aes(x=PC1,y=PC2, color = group)) +
  theme_classic(base_size = 20) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="Gene Expression", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs

# PCA: Heat stress treatment and heat stress day
deseq_PCA_hs_and_day <- plotPCA(vtd,intgroup=c("condition.stress.txt", ), returnData = TRUE) 
deseq_PCA_hs_and_day$group<- as.factor(deseq_PCA_hs_and_day$group)
deseq_PCA_hs_and_day$stress.day<- my.metadata$stress.day[match((deseq_PCA_hs_and_day$name), my.metadata$ID)]

# Plot PCA heat stress treatment
PCA_plot_hs_and_day <- ggplot(deseq_PCA_hs_and_day,aes(x=PC1,y=PC2, shape = group, color = stress.day)) +
  theme_classic(base_size = 20) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  # labs(title="Gene Expression", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs_and_day

# PCA: Heat stress treatment and acclimation treatment
deseq_PCA_hs_and_acclimation <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs_and_acclimation$group<- as.factor(deseq_PCA_hs_and_acclimation$group)
deseq_PCA_hs_and_acclimation$acclimation <-my.metadata$acclimation.group[match((deseq_PCA_hs_and_acclimation$name), my.metadata$ID)]

# Plot PCA heat stress treatment
PCA_plot_hs_and_acclimation <- ggplot(deseq_PCA_hs_and_acclimation,aes(x=PC1,y=PC2, shape = group, color = acclimation)) +
  theme_classic(base_size = 20) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  # labs(title="Gene Expression", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs_and_acclimation
