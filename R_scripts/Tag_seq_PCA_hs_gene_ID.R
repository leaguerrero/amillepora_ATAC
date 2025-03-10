## First round data exploration

# Download packages:
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# File Set-up
file.dir <- "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Gene_count_files/"
my.files <- grep(".txt", list.files(file.dir), value=TRUE)
my.metadata <- read.table("/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)
sum(my.metadata$acclimation.class == "acclimated")
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
                             condition = my.metadata[,2:10])

# Convert variables to factors
factorVars <- c("condition.acclimation.temp", "condition.acclimation.class", "condition.acclimation.txt",
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
    condition.acclimation.txt:condition.stress.txt)

# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


ddsHTSeq$condition.acclimation.txt <- relevel(ddsHTSeq$condition.acclimation.txt, "control")
ddsHTSeq$condition.acclimation.class <- relevel(ddsHTSeq$condition.acclimation.class, "control")

# Extract normalized counts matrix
# Estimate size factors
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)

# Estimate normalization factors
ddsHTSeq <- estimateDispersions(ddsHTSeq)
normalized.counts <- as.data.frame(counts(ddsHTSeq, normalized = TRUE))
normalized.counts$gene <- rownames(normalized.counts) 
normalized.counts$standard.deviation <- apply(normalized.counts[,1:58], 1, sd)


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
#####
# PCA: heat stress day
deseq_PCA_heat_stress_day <- plotPCA(vtd,intgroup=c( "condition.stress.day"), returnData = TRUE) 
deseq_PCA_heat_stress_day$condition.stress.day

# Plot PCA Heat Stress day
PCA_plot_hs_day <- ggplot(deseq_PCA_heat_stress_day,aes(x=PC1,y=PC2, color=group)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression colored by heat assay day") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs_day
# ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_hs_day.jpg", 
#        plot = PCA_plot_hs_day, 
#        width = 8, 
#        height = 6, 
#        dpi = 300,
#        device = "jpeg"
# )

# PCA: extraction batch
deseq_PCA_batch <- plotPCA(vtd,intgroup=c( "condition.ext.batch"), returnData = TRUE) 
deseq_PCA_batch$group<- as.factor(deseq_PCA_batch$group)

# Plot PCA ext batch
PCA_plot_ext_batch <- ggplot(deseq_PCA_batch,aes(x=PC1,y=PC2, color = group)) +
  theme_classic(base_size = 10) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="PCA of gene expression colored by extraction batch", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_ext_batch

# ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_extraction_batch.jpg", 
#        plot = PCA_plot_ext_batch, 
#        width = 8, 
#        height = 6, 
#        dpi = 300,
#        device = "jpeg"
# )
# PCA: Heat stress treatment
deseq_PCA_hs <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs$group<- as.factor(deseq_PCA_hs$group)

# Plot PCA heat stress treatment
PCA_plot_hs <- ggplot(deseq_PCA_hs,aes(x=PC1,y=PC2, color = condition.stress.txt)) +
  theme_classic(base_size = 20) +
  #theme(legend.position="none") + 
  scale_color_manual(values = c("#bfd7ea", "#ff5a5f"), 
                     name = "Treatment", 
                     labels = c("24.5°C", "36°C")) +
  labs(title="PCA of gene expression colored by heat stress treatment") +
  geom_point(size = 5, alpha = 0.8)

PCA_plot_hs
# ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_heat_vs_control.jpg", 
#        plot = PCA_plot_hs, 
#        width = 8, 
#        height = 6, 
#        dpi = 300,
#        device = "jpeg"
# )


patchwork <- (heat.response.boxplot + R.chan.intensity.bxp + glut.red.boxplot + PCA_plot_hs) 
patchwork + plot_annotation(
  title = expression(paste("Induced Heat Response of ", italic("A. millepora"))),
  tag_levels = list(c('A.', 'B.','C.', 'D.'))
)

# ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/Figure_2.jpg",
#        plot = patchwork + plot_annotation(
#          title = expression(paste("Induced Heat Response of ", italic("A. millepora"))),
#          tag_levels = list(c('A.', 'B.','C.', 'D.'))
#        ),
#        width = 12,
#        height = 10,
#        dpi = 300,
#        device = "jpeg"
# )


# PCA: Heat stress treatment and heat stress day
deseq_PCA_hs_and_day <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
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

# PCA: Heat stress treatment and acclimation class
deseq_PCA_hs_and_acclimation <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs_and_acclimation$group<- as.factor(deseq_PCA_hs_and_acclimation$group)
deseq_PCA_hs_and_acclimation$acclimation <-my.metadata$acclimation.class[match((deseq_PCA_hs_and_acclimation$name), my.metadata$ID)]

# Plot PCA heat stress treatment
PCA_plot_hs_and_acclimation <- ggplot(deseq_PCA_hs_and_acclimation,aes(x=PC1,y=PC2, shape = group, color = acclimation)) +
  theme_classic(base_size = 10) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="Gene Expression PCA", shape ="Heat stress treatment") +
  geom_point(size = 5, alpha = 0.7)
PCA_plot_hs_and_acclimation

ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_hs_and_acclimation.jpg", 
       plot = PCA_plot_hs_and_acclimation, 
       width = 8, 
       height = 6, 
       dpi = 300,
       device = "jpeg"
)

# Separation along PC2. Is this driven by genotype?
# PCA: Heat stress treatment and acclimation treatment
deseq_PCA_hs_and_genotype <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs_and_genotype$group<- as.factor(deseq_PCA_hs_and_genotype$group)
deseq_PCA_hs_and_genotype$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_hs_and_genotype$name), my.metadata$ID)])

# Plot PCA heat stress treatment
PCA_plot_hs_and_genotype <- ggplot(deseq_PCA_hs_and_genotype,aes(x=PC1,y=PC2, shape = group, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  labs(title="Gene Expression PCA", shape ="Heat stress treatment") +
  geom_point(size = 5) 
PCA_plot_hs_and_genotype

ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_hs_and_genotype.jpg", 
       plot = PCA_plot_hs_and_genotype, 
       width = 8, 
       height = 6, 
       dpi = 300,
       device = "jpeg"
)


## Look at heat stress and acclimation group (instead of acclimated and non acclimated)
# Plot PCA heat stress treatment
deseq_PCA_hs_and_acclimation_treatment <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE) 
deseq_PCA_hs_and_acclimation_treatment$group<- as.factor(deseq_PCA_hs_and_acclimation_treatment$group)
deseq_PCA_hs_and_acclimation_treatment$acclimation.treatment <-my.metadata$acclimation.txt[match((deseq_PCA_hs_and_acclimation_treatment$name), my.metadata$ID)]


PCA_plot_hs_and_acclimation_treatment <- ggplot(deseq_PCA_hs_and_acclimation_treatment,aes(x=PC1,y=PC2, shape = group, color = acclimation.treatment)) +
  theme_classic(base_size = 14) +
  #theme(legend.position="none") + 
  scale_color_manual(values = c("#F6D97C", "#F3C995", "#ECAAC6", "#E68BF8")) +
  labs(title="Gene Expression PCA", shape ="Heat stress treatment") +
  geom_point(size = 7)
PCA_plot_hs_and_acclimation_treatment

ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_hs_and_acclimation.jpg", 
       plot = PCA_plot_hs_and_acclimation, 
       width = 8, 
       height = 6, 
       dpi = 300,
       device = "jpeg"
)

# PC1 vs PC3
#########################
deseq_PCA_hs_and_acclimation_treatment <- plotPCA(vtd,intgroup=c("condition.stress.txt"), returnData = TRUE, pcsToUse = c(3,4)) 
deseq_PCA_hs_and_acclimation_treatment$group<- as.factor(deseq_PCA_hs_and_acclimation_treatment$group)
deseq_PCA_hs_and_acclimation_treatment$acclimation.treatment <-my.metadata$acclimation.txt[match((deseq_PCA_hs_and_acclimation_treatment$name), my.metadata$ID)]


PCA_plot_hs_and_acclimation_treatment <- ggplot(deseq_PCA_hs_and_acclimation_treatment,aes(x=PC3,y=PC4, shape = group, color = acclimation.treatment)) +
  theme_classic(base_size = 14) +
  #theme(legend.position="none") + 
  scale_color_manual(values = c("#F6D97C", "#F3C995", "#ECAAC6", "#E68BF8")) +
  labs(title="Gene Expression PCA", shape ="Heat stress treatment") +
  geom_point(size = 7)
PCA_plot_hs_and_acclimation_treatment

# PCA: heat stress day
deseq_PCA_heat_stress_day <- plotPCA(vtd,intgroup=c( "condition.stress.day"), returnData = TRUE, pcsToUse = c(3,4)) 
deseq_PCA_heat_stress_day$condition.stress.day

# Plot PCA Heat Stress day
PCA_plot_hs_day <- ggplot(deseq_PCA_heat_stress_day,aes(x=PC3,y=PC4, color=group)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression colored by heat assay day") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_hs_day

## PC 1 vs PC 3 batch
deseq_PCA_batch <- plotPCA(vtd,intgroup=c( "condition.ext.batch"), returnData = TRUE, pcsToUse = c(1,3)) 
deseq_PCA_batch$group<- as.factor(deseq_PCA_batch$group)


PCA_plot_ext_batch <- ggplot(deseq_PCA_batch,aes(x=PC1,y=PC3, color = group)) +
  theme_classic(base_size = 10) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="PCA of gene expression colored by extraction batch", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_ext_batch

#####
# Identifying heat stress genes:
## Look at only control sample heat stress genes since acclimation can dampen the ge response
non.acclimated.files.list <- filter(my.metadata, acclimation.class == "control")[,11]

# Filter metadata and create sample table
fil.metadata <- filter(my.metadata, acclimation.class == "control")
non.acc.sampleNames <- fil.metadata$ID

non.acclimated.sampleTable <- data.frame(
  sampleName = non.acc.sampleNames,
  fileName = non.acclimated.files.list,
  condition = fil.metadata[,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp", 
                          "condition.stress.txt",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  non.acclimated.sampleTable[[column]] <- factor(non.acclimated.sampleTable[[column]])
}


# Create DESeqDataSet
non_acclimated_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = non.acclimated.sampleTable,
  directory = file.dir,
  design = ~ condition.stress.txt
)

# Filter low count rows
nonacc.keep <- rowSums(counts(non_acclimated_ddsHTSeq)) >= 10
non_acclimated_ddsHTSeq <- non_acclimated_ddsHTSeq[nonacc.keep,]

# Differential Expression Analysis
non.acc.dds <- DESeq(non_acclimated_ddsHTSeq)
resultsNames(non.acc.dds)
levels(non.acc.dds$condition.stress.txt)

# heat stress vs control
## Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq 
## run without LFC shrinkage. For consistency with results, the column name 
## lfcSE is used here although what is returned is a posterior SD. 
names(non.acc.dds@colData) # condition.stress.txt
hs.results <- lfcShrink(non.acc.dds,
                        contrast = c("condition.stress.txt", "36", "24.5"),
                        type = "ashr")

heat.response.results <- hs.results[complete.cases(hs.results),]

nonacc.HRGs <- subset(heat.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

nonacc.HRGs.tb <- nonacc.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
dim(nonacc.HRGs.tb)
# 2212    6

# write.csv(nonacc.HRGs.tb, file = "~/Lab Notebook/Chapter2/Data_analysis/TagSeq/Data_out/
", row.names = FALSE)
#####

## Look at non-heat stressed samples PCA
non.hs.files.list <- filter(my.metadata, stress.group == "control")[,11]

# Filter metadata and create sample table
no.hs.metadata <- filter(my.metadata, stress.group == "control")
no.hs.sampleNames <- no.hs.metadata$ID

no.hs.sampleTable <- data.frame(
  sampleName = no.hs.sampleNames,
  fileName = non.hs.files.list,
  condition = no.hs.metadata [,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.txt",
                          "condition.acclimation.class",
                         "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  no.hs.sampleTable[[column]] <- factor(no.hs.sampleTable[[column]])
}


# Create DESeqDataSet
no_hs_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = no.hs.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.txt + condition.genotype 
)

no_hs_ddsHTSeq$condition.acclimation.txt <- relevel(no_hs_ddsHTSeq$condition.acclimation.txt, "control")
# Filter low count rows
no_hs_keep <- rowSums(counts(no_hs_ddsHTSeq)) >= 10
no_hs_ddsHTSeq<- no_hs_ddsHTSeq[no_hs_keep,]

# Differential Expression Analysis
no.hs.dds <- DESeq(no_hs_ddsHTSeq)
resultsNames(no.hs.dds)
levels(no.hs.dds$condition.acclimation.class)

# No heat stress PCA
vtd.no.hs <- vst(no.hs.dds)

# PCA: GE of control samples colored by acclimation day
deseq_PCA_acclimation_day <- plotPCA(vtd.no.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
deseq_PCA_acclimation_day$hs.day <-as.factor(my.metadata$stress.day[match((deseq_PCA_acclimation_day$name), my.metadata$ID)])

#deseq_PCA_acclimation_day$condition.acclimation.class

# Plot PCA Heat Stress day


PCA_plot_hs_day <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color=hs.day)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression (control samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_hs_day

# PCA of acclimated controls colored by genotype:
deseq_PCA_acclimation_day <- plotPCA(vtd.no.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
deseq_PCA_acclimation_day$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day$name), my.metadata$ID)])

# Plot PCA heat stress treatment
PCA_plot_acc_and_genotype <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
PCA_plot_acc_and_genotype

PCA_plot_acc_day + PCA_plot_acc_and_genotype


# PCA 2 & 3 of acclimated controls colored by genotype:
deseq_PCA_acclimation_day_2_3 <- plotPCA(vtd.no.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE, pcsToUse = c(2,3)) 
deseq_PCA_acclimation_day_2_3$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day_2_3$name), my.metadata$ID)])
deseq_PCA_acclimation_day_2_3$batch <-as.factor(my.metadata$ext.batch[match((deseq_PCA_acclimation_day_2_3$name), my.metadata$ID)])
deseq_PCA_acclimation_day_2_3$hs.day <-as.factor(my.metadata$stress.day[match((deseq_PCA_acclimation_day_2_3$name), my.metadata$ID)])

# Plot PCA 
PCA_plot_acc_and_genotype_2_3_geno <- ggplot(deseq_PCA_acclimation_day_2_3,aes(x=PC2,y=PC3, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  geom_point(size = 5) +
  theme_pubr()
PCA_plot_acc_and_genotype_2_3_geno

PCA_plot_acc_and_genotype_2_3 <- ggplot(deseq_PCA_acclimation_day_2_3,aes(x=PC2,y=PC3, color = hs.day)) +
  theme_classic(base_size = 10) +
  geom_point(size = 5) +
  theme_pubr()
PCA_plot_acc_and_genotype_2_3

PCA_plot_acc_and_genotype_2_3 + PCA_plot_acc_and_genotype_2_3_geno

# Batch effects?
# PCA: extraction batch
no_hs_deseq_PCA_batch <- plotPCA(vtd.no.hs,intgroup=c( "condition.ext.batch"), returnData = TRUE, pcsToUse = c(2,3)) 
no_hs_deseq_PCA_batch$group<- as.factor(no_hs_deseq_PCA_batch$group)

# Plot PCA ext batch
PCA_plot_ext_batch_no_hs <- ggplot(no_hs_deseq_PCA_batch,aes(x=PC2,y=PC3, color = group)) +
  theme_classic(base_size = 10) +
  #theme(legend.position="none") + 
  #scale_color_manual(values = c("#66BBBB", "#D12E8883")) +
  labs(title="PCA of gene expression colored by extraction batch", shape ="Ext Batch") +
  geom_point(size = 5, alpha = 0.7)

PCA_plot_ext_batch_no_hs



# Look at only HS samples

hs.files.list <- filter(my.metadata, stress.group == "heat")[,11]

# Filter metadata and create sample table
hs.metadata <- filter(my.metadata, stress.group == "heat")
hs.sampleNames <- hs.metadata$ID

hs.sampleTable <- data.frame(
  sampleName = hs.sampleNames,
  fileName = hs.files.list,
  condition = hs.metadata [,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.txt",
                          "condition.acclimation.class",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  hs.sampleTable[[column]] <- factor(hs.sampleTable[[column]])
}


# Create DESeqDataSet
hs_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = hs.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.txt + condition.genotype
)

hs_ddsHTSeq$condition.acclimation.txt <- relevel(hs_ddsHTSeq$condition.acclimation.txt, "control")
# Filter low count rows
hs_keep <- rowSums(counts(hs_ddsHTSeq)) >= 10
hs_ddsHTSeq<- hs_ddsHTSeq[hs_keep,]

# Differential Expression Analysis
hs.dds <- DESeq(hs_ddsHTSeq)
resultsNames(hs.dds)
levels(hs.dds$condition.acclimation.class)

# No heat stress PCA
vtd.hs <- vst(hs.dds)

# PCA: GE of control samples colored by acclimation day
deseq_PCA_acclimation_day_hs <- plotPCA(vtd.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
#deseq_PCA_acclimation_day$condition.acclimation.class

# Plot PCA Heat Stress day
PCA_plot_acc_day_hs <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color=group)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression colored by acclimation treatment") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc_day_hs 

deseq_PCA_acclimation_day_hs <- plotPCA(vtd.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
deseq_PCA_acclimation_day_hs$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day_hs$name), my.metadata$ID)])
deseq_PCA_acclimation_day_hs$batch <-as.factor(my.metadata$ext.batch[match((deseq_PCA_acclimation_day_hs$name), my.metadata$ID)])


deseq_PCA_acclimation_day_hs <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
deseq_PCA_acclimation_day_hs

PCA_plot_acc_day_hs + deseq_PCA_acclimation_day_hs


deseq_PCA_acclimation_day_hs_2_3 <- plotPCA(vtd.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE, pcsToUse = c(2,3)) 
deseq_PCA_acclimation_day_hs_2_3$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day_hs_2_3$name), my.metadata$ID)])
deseq_PCA_acclimation_day_hs_2_3$batch <-as.factor(my.metadata$ext.batch[match((deseq_PCA_acclimation_day_hs_2_3$name), my.metadata$ID)])

deseq_PCA_acclimation_day_hs_2_3_acc <- ggplot(deseq_PCA_acclimation_day_hs_2_3,aes(x=PC2,y=PC3, color = group)) +
  theme_classic(base_size = 10) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
deseq_PCA_acclimation_day_hs_2_3_acc

deseq_PCA_acclimation_day_hs_2_3 <- ggplot(deseq_PCA_acclimation_day_hs_2_3,aes(x=PC2,y=PC3, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
deseq_PCA_acclimation_day_hs_2_3

PCA_plot_acc_day_hs + deseq_PCA_acclimation_day_hs




## ID transcriptionally modulated genes on D7
D07.files.list <- filter(my.metadata,stress.day == "D07")[,11]

# Filter metadata and create sample table
D07.metadata <- filter(my.metadata,stress.day == "D07")
D07.sampleNames <- D07.metadata$ID

D07.sampleTable <- data.frame(
  sampleName = D07.sampleNames,
  fileName = D07.files.list,
  condition = D07.metadata[,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.class",
                          "condition.stress.group",
                          "condition.stress.txt",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  D07.sampleTable[[column]] <- factor(D07.sampleTable[[column]])
}



# Create DESeqDataSet Object
D7ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = D07.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.class +
    condition.stress.group +
    condition.stress.group:condition.acclimation.class)
D7ddsHTSeq$condition.acclimation.class <- relevel(D7ddsHTSeq$condition.acclimation.class, "control")
# Pre-filter
keep <- rowSums(counts(D7ddsHTSeq)) >= 10
D7ddsHTSeq <- D7ddsHTSeq[keep,]

# Differential Expression Analysis
D7.dds <- DESeq(D7ddsHTSeq) 
resultsNames(D7.dds)


D7.coeff.df <- as.data.frame(coef(D7.dds), colnames = TRUE)
D7.wald.stat.df <- as.data.frame(cbind(row.names(D7.coeff.df), as.numeric(D7.coeff.df$condition.acclimation.classacclimated.condition.stress.groupheat)), colnames = TRUE)
colnames(D7.wald.stat.df) <- c("gene", "coefficient")

D7.significant.results <- results(D7.dds, name = "condition.acclimation.classacclimated.condition.stress.groupheat")
D7.significant.results <- D7.significant.results[complete.cases(D7.significant.results),]
D7.significant.results <- D7.significant.results[D7.significant.results$padj < 0.1,]
D7.significant.results.df <- D7.significant.results %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") 

# significant.interaction.df <- dplyr::left_join(significant.results.df, 
#                                                wald.stat.df, 
#                                                by = "gene") %>% 
#   mutate(coefficient = as.numeric(coefficient))



# Dampened Genes
neg.wald.stat <- nrow(significant.interaction.df[significant.interaction.df$stat < 0,])
dampened.genes <- significant.interaction.df[significant.interaction.df$stat < 0,]
# write.csv(dampened.genes, file = "~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.csv", row.names = FALSE)

# Amplified Genes
pos.wald.stat <- nrow(significant.interaction.df[significant.interaction.df$stat > 0,])
amp.genes <- significant.interaction.df[significant.interaction.df$stat > 0,]
# write.csv(amp.genes, file = "~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.csv", row.names = FALSE)


## ID transcriptionally modulated genes on D14
D14.files.list <- filter(my.metadata,stress.day == "D14")[,11]

# Filter metadata and create sample table
D14.metadata <- filter(my.metadata,stress.day == "D14")
D14.sampleNames <- D14.metadata$ID

D14.sampleTable <- data.frame(
  sampleName = D14.sampleNames,
  fileName = D14.files.list,
  condition = D14.metadata[,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.class",
                          "condition.stress.group",
                          "condition.stress.txt",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  D14.sampleTable[[column]] <- factor(D14.sampleTable[[column]])
}



# Create DESeqDataSet Object
D14ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = D14.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.class +
    condition.stress.group +
    condition.stress.group:condition.acclimation.class)
D14ddsHTSeq$condition.acclimation.class <- relevel(D14ddsHTSeq$condition.acclimation.class, "control")
# Pre-filter
rm(keep)
keep <- rowSums(counts(D14ddsHTSeq)) >= 10
D14ddsHTSeq <- D14ddsHTSeq[keep,]

# Differential Expression Analysis
D14.dds <- DESeq(D14ddsHTSeq) 
resultsNames(D14.dds)


D14.coeff.df <- as.data.frame(coef(D14.dds), colnames = TRUE)
D14.wald.stat.df <- as.data.frame(cbind(row.names(D14.coeff.df), as.numeric(D14.coeff.df$condition.acclimation.classacclimated.condition.stress.groupheat)), colnames = TRUE)
colnames(D14.wald.stat.df) <- c("gene", "coefficient")

D14.significant.results <- results(D14.dds, name = "condition.acclimation.classacclimated.condition.stress.groupheat")
D14.significant.results <- D14.significant.results[complete.cases(D14.significant.results),]
D14.significant.results <- D14.significant.results[D14.significant.results$padj < 0.1,]
D14.significant.results.df <- D14.significant.results %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") 

## Heat stress genes D0 only:
D0.files.list <- filter(my.metadata, stress.day == "D00")[,11]

# Filter metadata and create sample table
D0.metadata <- filter(my.metadata, stress.day == "D00")
D0.sampleNames <- D0.metadata$ID

D0.sampleTable <- data.frame(
  sampleName = D0.sampleNames,
  fileName = D0.files.list,
  condition = D0.metadata[,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.stress.txt",
                          "condition.stress.group",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
 D0.sampleTable[[column]] <- factor(D0.sampleTable[[column]])
}


# Create DESeqDataSet
D0_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = D0.sampleTable,
  directory = file.dir,
  design = ~ condition.stress.group
)

# Filter low count rows
D0.keep <- rowSums(counts(D0_ddsHTSeq)) >= 10
D0_ddsHTSeq <- D0_ddsHTSeq[D0.keep,]

# counts matrix
D0.normalized.counts <- counts(D0_ddsHTSeq, normalized = TRUE)

# Differential Expression Analysis
D0.hs.dds <- DESeq(D0_ddsHTSeq)
resultsNames(D0.hs.dds)
levels(D0.hs.dds$condition.stress.group)

# heat stress vs control
## Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq 
## run without LFC shrinkage. For consistency with results, the column name 
## lfcSE is used here although what is returned is a posterior SD. 
names(D0.hs.dds@colData) # condition.stress.txt
D0.hs.results <- lfcShrink(D0.hs.dds,
                        contrast = c("condition.stress.group", "heat", "control"),
                        type = "ashr")

D0.heat.response.results <- D0.hs.results[complete.cases(D0.hs.results),]

D0.HRGs <- subset(D0.heat.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

D0.HRGs.tb <- D0.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
dim(D0.HRGs.tb)
# 1054    6

# write.csv(D0.HRGs.tb, file = '~/Lab Notebook/Chapter2/Data_analysis/TagSeq/Data_out/D0.heat.response.genes.csv', row.names = FALSE)

## PCA From heat stress genes only
## Control Samples
## Look at non-heat stressed samples PCA
non.hs.files.list <- filter(my.metadata, stress.group == "control")[,11]

# Filter metadata and create sample table
no.hs.metadata <- filter(my.metadata, stress.group == "control")
no.hs.sampleNames <- no.hs.metadata$ID

no.hs.sampleTable <- data.frame(
  sampleName = no.hs.sampleNames,
  fileName = non.hs.files.list,
  condition = no.hs.metadata [,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.txt",
                          "condition.acclimation.class",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  no.hs.sampleTable[[column]] <- factor(no.hs.sampleTable[[column]])
}


# Create DESeqDataSet

no_hs_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = no.hs.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.txt + condition.genotype 
)

no_hs_ddsHTSeq$condition.acclimation.txt <- relevel(no_hs_ddsHTSeq$condition.acclimation.txt, "control")
# Filter low count rows
no_hs_keep <- rowSums(counts(no_hs_ddsHTSeq)) >= 10
no_hs_ddsHTSeq<- no_hs_ddsHTSeq[no_hs_keep,]

# Differential Expression Analysis
no.hs.dds <- DESeq(no_hs_ddsHTSeq)
resultsNames(no.hs.dds)
levels(no.hs.dds$condition.acclimation.class)

# No heat stress PCA
D0.heat.response.genes <- D0.HRGs.tb$gene
vtd.no.hs <- vst(no.hs.dds)
filtered.vtd.no.hs <- vtd.no.hs[rownames(rowData(vtd.no.hs))%in% D0.heat.response.genes, ]
# PCA: GE of control samples colored by acclimation day
deseq_PCA_acclimation_day <- plotPCA(filtered.vtd.no.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
deseq_PCA_acclimation_day$hs.day <-as.factor(my.metadata$stress.day[match((deseq_PCA_acclimation_day$name), my.metadata$ID)])
deseq_PCA_acclimation_day$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day$name), my.metadata$ID)])
deseq_PCA_acclimation_day$acclimation.temp <-as.factor(my.metadata$acclimation.temp[match((deseq_PCA_acclimation_day$name), my.metadata$ID)])
#deseq_PCA_acclimation_day$condition.acclimation.class

# Plot PCA Heat Stress day

PCA_plot_hs_day <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color=hs.day, shape = acclimation.temp)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of heat response gene expression (control samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_hs_day

PCA_plot_acc <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color=group)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of heat response gene expression (control samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc


PCA_plot_acc_hs <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color=group, shape = hs.day)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of heat response gene expression (control samples)") +
  scale_shape_manual(values = c(17, 19, 15,0)) +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc_hs

deseq_PCA_genotype <- ggplot(deseq_PCA_acclimation_day,aes(x=PC1,y=PC2, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
deseq_PCA_genotype

PCA_plot_hs_day + PCA_plot_acc 



### Heat stress samples looking at how samples cluster based on heat response genes
hs.files.list <- filter(my.metadata, stress.group == "heat")[,11]

# Filter metadata and create sample table
hs.metadata <- filter(my.metadata, stress.group == "heat")
hs.sampleNames <- hs.metadata$ID

hs.sampleTable <- data.frame(
  sampleName = hs.sampleNames,
  fileName = hs.files.list,
  condition = hs.metadata [,2:10]
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.acclimation.temp",
                          "condition.acclimation.txt",
                          "condition.acclimation.class",
                          "condition.ext.batch",
                          "condition.genotype")

for (column in columns_to_factorize) {
  hs.sampleTable[[column]] <- factor(hs.sampleTable[[column]])
}


# Create DESeqDataSet
hs_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = hs.sampleTable,
  directory = file.dir,
  design = ~ condition.acclimation.txt + condition.genotype
)

hs_ddsHTSeq$condition.acclimation.txt <- relevel(hs_ddsHTSeq$condition.acclimation.txt, "control")
# Filter low count rows
hs_keep <- rowSums(counts(hs_ddsHTSeq)) >= 10
hs_ddsHTSeq<- hs_ddsHTSeq[hs_keep,]

# Differential Expression Analysis
hs.dds <- DESeq(hs_ddsHTSeq)
resultsNames(hs.dds)
levels(hs.dds$condition.acclimation.class)

# No heat stress PCA
vtd.hs <- vst(hs.dds)
filtered.vtd.hs <- vtd.hs[rownames(rowData(vtd.hs))%in% D0.heat.response.genes, ]

# PCA: GE of control samples colored by acclimation day
deseq_PCA_acclimation_day_hs <- plotPCA(filtered.vtd.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE) 
deseq_PCA_acclimation_day_hs$hs.day <-as.factor(my.metadata$stress.day[match((deseq_PCA_acclimation_day_hs$name), my.metadata$ID)])
deseq_PCA_acclimation_day_hs$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day_hs$name), my.metadata$ID)])

#deseq_PCA_acclimation_day$condition.acclimation.class

# Plot PCA Heat Stress day

PCA_plot_hs_day <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color=hs.day)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression (heated samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_hs_day

PCA_plot_acc <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color=group)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression (heated samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc


deseq_PCA_genotype <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color = genotype)) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74","#EFB366", "#EE8959", "#E76F51")) +
  #labs(title="Gene Expression PCA") +
  geom_point(size = 5)+
  theme_pubr()
deseq_PCA_genotype


PCA_plot_acc_hs_hrg <- ggplot(deseq_PCA_acclimation_day_hs,aes(x=PC1,y=PC2, color=group, shape = hs.day)) +
  theme_classic(base_size = 10) +
  #labs(title="PCA of heat response gene expression (control samples)") +
  scale_shape_manual(values = c(17, 19, 15,0)) +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc_hs_hrg


# PCA: GE of control samples colored by acclimation day
deseq_PCA_acclimation_day_hs_2_3 <- plotPCA(filtered.vtd.hs,intgroup=c( "condition.acclimation.txt"), returnData = TRUE, pcsToUse = c(3,4)) 
deseq_PCA_acclimation_day_hs_2_3$hs.day <-as.factor(my.metadata$stress.day[match((deseq_PCA_acclimation_day_hs_2_3$name), my.metadata$ID)])
deseq_PCA_acclimation_day_hs_2_3$genotype <-as.factor(my.metadata$genotype[match((deseq_PCA_acclimation_day_hs_2_3$name), my.metadata$ID)])
deseq_PCA_acclimation_day_hs_2_3$acclimation.temp <-as.factor(my.metadata$acclimation.temp[match((deseq_PCA_acclimation_day_hs_2_3$name), my.metadata$ID)])

#deseq_PCA_acclimation_day$condition.acclimation.class

# Plot PCA Heat Stress day

PCA_plot_hs_day <- ggplot(deseq_PCA_acclimation_day_hs_2_3,aes(x=PC3,y=PC4, color=hs.day, shape= acclimation.temp)) +
  theme_classic(base_size = 10) +
  labs(title="PCA of gene expression (heated samples)") +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_hs_day


PCA_plot_acc_hs_hrg <- ggplot(deseq_PCA_acclimation_day_hs_2_3,aes(x=PC3,y=PC4, color=group, shape = hs.day)) +
  theme_classic(base_size = 10) +
  #labs(title="PCA of heat response gene expression (control samples)") +
  scale_shape_manual(values = c(17, 19, 15,0)) +
  geom_point(size = 5) +
  theme_pubr()

PCA_plot_acc_hs_hrg







