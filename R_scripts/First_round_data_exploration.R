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
  theme_classic(base_size = 10) +
  #theme(legend.position="none") + 
  scale_color_manual(values = c("#a7b7c4", "#F6D97C"), 
                     name = "Heat Stress Treatment") + #, 
                   #  labels = c("24.5°C", "36°C")) +
  labs(title="PCA of gene expression colored by heat stress treatment") +
  geom_point(size = 5)

PCA_plot_hs
# ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_heat_vs_control.jpg", 
#        plot = PCA_plot_hs, 
#        width = 8, 
#        height = 6, 
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

hs.results <- lfcShrink(non.acc.dds,
                        contrast = "condition.stress.txt_36_vs_24.5", "36", "24.5",
                        type = "ashr")
D0.29.response.results <- D0.29.response.results[complete.cases(D0.29.response.results),]

D0.29C.HRGs <- subset(D0.29.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

D0.29C.HRGs.tb <- D0.29C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Day 0 31
D0.31.response.results <- lfcShrink(group.dds,
                                    contrast = c("group", "stress031", "control031"),
                                    type = "ashr")
D0.31.response.results <- D0.31.response.results[complete.cases(D0.31.response.results),]

D0.31C.HRGs <- subset(D0.31.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

D0.31C.HRGs.tb <- D0.31C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Day 11 29
D11.29.response.results <- lfcShrink(group.dds,
                                     contrast = c("group", "stress1129", "control1129"),
                                     type = "ashr")
D11.29.response.results <- D11.29.response.results[complete.cases(D11.29.response.results),]

D11.29C.HRGs <- subset(D11.29.response.results,
                       padj < 0.01 &
                         abs(log2FoldChange) >= 2)

D11.29C.HRGs.tb <- D11.29C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
dim(D11.29C.HRGs.tb)

nonaccl.HRG <- union(union(D0.31C.HRGs.tb$gene, 
                           D11.29C.HRGs.tb$gene), 
                     D0.29C.HRGs.tb$gene)
nonaccl.HRG.gene.list <- as.data.frame(nonaccl.HRG)
colnames(nonaccl.HRG.gene.list) <- 'gene'
# write.csv(nonaccl.HRG.gene.list, file = '~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/nonaccl.HRG.gene.list.csv', row.names = FALSE)


core.HRGs <- inner_join(
  inner_join(D0.31C.HRGs.tb, 
             D11.29C.HRGs.tb, 
             by = "gene"), 
  D0.29C.HRGs.tb, by = "gene")


core.HRGs.gene.list <- as.data.frame(core.HRGs$gene)
colnames(core.HRGs.gene.list) <- 'gene'
# write.csv(core.HRGs.gene.list, file = '~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.HRG.list.csv', row.names = FALSE)



