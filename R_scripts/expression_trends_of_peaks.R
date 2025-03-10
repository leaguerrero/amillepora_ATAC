# Chapter 2: Peaks over constitutively expressed genes

# Load libraries
library(GenomicRanges)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(ggVennDiagram)
library(beeswarm)
library(ggdist)
library(gghalves)
library(ggpubr)

# Load in data sets
# Transcripts
gff.file <- "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amil.all.maker.noseq.gff"
txdb <- makeTxDbFromGFF(gff.file, format = "gff3") 
a.millepora.txs <- transcripts(txdb)

# Files of called peaks (MACS2 Output)
file_paths <- c("~/Lab Notebook/Chapter2/Data_analysis/ATAC/Peak_files/AMIL_GENO_01_S1_summits.bed",
                "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Peak_files/AMIL_GENO_02_S2_summits.bed",
                "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Peak_files/AMIL_GENO_03_S3_summits.bed",
                "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Peak_files/AMIL_GENO_04_S4_summits.bed")

peak.gr.list <- list()

for (file_path in file_paths) {
  peaks <- readPeakFile(file_path)
  peak.gr.list[[file_path]] <- peaks
}

# Peak Annotation
peakAnno.1 <- annotatePeak(peak.gr.list[[1]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

peakAnno.2 <- annotatePeak(peak.gr.list[[2]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

peakAnno.3 <- annotatePeak(peak.gr.list[[3]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

peakAnno.4 <- annotatePeak(peak.gr.list[[4]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

# Annotated peaks data frames
peaks.df.1 <- as.data.frame(peakAnno.1@anno@elementMetadata@listData)
genomic.annotation.detail.df.1 <- peakAnno.1@detailGenomicAnnotation
peaks.df.1 <- cbind(peaks.df.1, genomic.annotation.detail.df.1)

peaks.df.2 <- as.data.frame(peakAnno.2@anno@elementMetadata@listData)
genomic.annotation.detail.df.2 <- peakAnno.2@detailGenomicAnnotation
peaks.df.2 <- cbind(peaks.df.2, genomic.annotation.detail.df.2)

peaks.df.3 <- as.data.frame(peakAnno.3@anno@elementMetadata@listData)
genomic.annotation.detail.df.3 <- peakAnno.3@detailGenomicAnnotation
peaks.df.3 <- cbind(peaks.df.3, genomic.annotation.detail.df.3)

peaks.df.4 <- as.data.frame(peakAnno.4@anno@elementMetadata@listData)
genomic.annotation.detail.df.4 <- peakAnno.4@detailGenomicAnnotation
peaks.df.4 <- cbind(peaks.df.4, genomic.annotation.detail.df.4)

# Isolate Promoter (<=1kb) separately across genotypes

geno.1.promoter <- filter(peaks.df.1, annotation == "Promoter (<=1kb)") # 3453

geno.2.promoter <- filter(peaks.df.2, annotation == "Promoter (<=1kb)") # 10540 

geno.3.promoter <- filter(peaks.df.3, annotation == "Promoter (<=1kb)") # 1274

geno.4.promoter <- filter(peaks.df.4, annotation == "Promoter (<=1kb)") # 3080

# Venn Diagram 

summary(geno.1.promoter$V5)
geno.1.first.quartile <- quantile(geno.1.promoter$V5, probs = 0.25)
geno.1.promoter <- filter(geno.1.promoter, V5 >= geno.1.first.quartile) #2616 promoters with peak score greater than the first quartile
genes.geno.1 <- data.frame(gene = geno.1.promoter$geneId)
genes.geno.1$geno<- as.factor("1")



summary(geno.2.promoter$V5)
geno.2.first.quartile <- quantile(geno.2.promoter$V5, probs = 0.25)
geno.2.promoter <- filter(geno.2.promoter, V5 >= geno.2.first.quartile) #8090 promoters with peak score greater than the first quartile
genes.geno.2 <- data.frame(gene = geno.2.promoter$geneId)
genes.geno.2$geno<- as.factor("2")



summary(geno.3.promoter$V5)
geno.3.first.quartile <- quantile(geno.3.promoter$V5, probs = 0.25)
geno.3.promoter <- filter(geno.3.promoter, V5 >= geno.3.first.quartile) #955 promoters with peak score greater than the first quartile
genes.geno.3 <- data.frame(gene = geno.3.promoter$geneId)
genes.geno.3$geno<- as.factor("3")

summary(geno.4.promoter$V5)
geno.4.first.quartile <- quantile(geno.4.promoter$V5, probs = 0.25)
geno.4.promoter <- filter(geno.4.promoter, V5 >= geno.4.first.quartile) #2341 promoters with peak score greater than the first quartile
genes.geno.4 <- data.frame(gene = geno.4.promoter$geneId)
genes.geno.4$geno<- as.factor("4")

x <- list(
  Geno_1 = genes.geno.1$gene, 
  Geno_2 = genes.geno.2$gene, 
  Geno_3 = genes.geno.3$gene,
  Geno_4 = genes.geno.4$gene
)
ggVennDiagram(x)

# Select genes that appear in at all 4 genos:
peaks.intersection.all.genos <- Reduce(intersect, x)

# Visualize the accessible regions and the expression in the control samples 
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
non.acc.dds <- DESeq(non_acclimated_ddsHTSeq)
non.acc.expression.df <- as.data.frame(results(non.acc.dds)) # Note:log2 fold change (MLE): condition.stress.txt 36 vs 24.5 
non.acc.expression.df <- non.acc.expression.df[non.acc.expression.df$padj < 0.05, ]

# Add peak data to non.acc.expression.df:
## 1st: TRUE or FALSE for if peaks are shared across all genotypes
in.peaks.all.genos <- rownames(non.acc.expression.df) %in% peaks.intersection.all.genos 

# Ground truth calls
true_genes <- non.acc.expression.df[in.peaks.all.genos, ]

# Assign TRUE/FALSE values based on presence in peaks.intersection.all.genos
non.acc.expression.df$all.genos <- in.peaks.all.genos

# Visualize
non.acc.expression.df$all.genos <- factor(non.acc.expression.df$all.genos)
filtered.non.acc.expression.df <- non.acc.expression.df[non.acc.expression.df$baseMean <= 20000,]
non.acc.expression.df$log10bM <- log10(non.acc.expression.df$baseMean)

if (any(non.acc.expression.df$baseMean <= 0)) {
  warning("Zero or negative values detected in baseMean column.")
}

# Filter out zero or negative values from baseMean
non.acc.expression.df <- non.acc.expression.df[non.acc.expression.df$baseMean > 0, ]

# Perform logarithm transformation
non.acc.expression.df$log10bM <- log10(non.acc.expression.df$baseMean)
non.acc.expression.df <- non.acc.expression.df[!is.na(non.acc.expression.df$log10bM), ]

# Check for NA values in log10bM
if (any(is.na(non.acc.expression.df$log10bM))) {
  warning("NA values detected in log10bM column.")
}


# Log10 Base mean expression of genes with accessible promoters
ggplot(non.acc.expression.df, aes(x = all.genos, 
                                  y = log10bM, 
                                  fill = all.genos)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#b5e2fa",
                               "#ed6a5a")) +
  labs(title = "Log10 base mean expression of genes with\naccessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 

# log2FC to heat stress of genes with accessible promoters
ggplot(non.acc.expression.df, aes(x = all.genos, 
                                  y = log2FoldChange, 
                                  fill = all.genos)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#b5e2fa",
                               "#ed6a5a")) +
  labs(title = "Gene log2FC in response to heat stress with\naccessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 



# lfcSE of accessible promoters
ggplot(non.acc.expression.df, aes(x = all.genos, 
                                  y = lfcSE, 
                                  fill = all.genos)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#b5e2fa",
                               "#ed6a5a")) +
  labs(title = "Gene expression variation with\naccessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 



# Peaks in all genos is too conservative to compare against genes with no peaks
## Look at peaks found in at least one genotype
## 2nd: TRUE or FALSE for if peak is found in at least one genotype
union.peaks <- Reduce(union, x) # 9859

in.peaks.one.geno <- rownames(non.acc.expression.df) %in% union.peaks #  genes overlap

# Ground truth calls
true_genes <- non.acc.expression.df[in.peaks.one.geno , ]

# Assign TRUE/FALSE values based on presence in peaks.intersection.all.genos
non.acc.expression.df$one.geno <- in.peaks.one.geno

# Visualize
non.acc.expression.df$one.geno <- factor(non.acc.expression.df$one.geno)

# Log10 Base mean expression of genes with accessible promoters
ggplot(non.acc.expression.df, aes(x = one.geno, 
                                  y = log10bM, 
                                  fill = one.geno)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#003566",
                               "#ada7c9")) +
  labs(title = "Expression of genes with accessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 

# log2FC to heat stress of genes with accessible promoters

ggplot(non.acc.expression.df, aes(x = one.geno, 
                                  y = log2FoldChange, 
                                  fill = one.geno)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#003566",
                               "#ada7c9")) +
  labs(title = "Heat stress l2fc in genes with accessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 



# lfcSE of accessible promoters
ggplot(non.acc.expression.df, aes(x = one.geno, 
                                  y = lfcSE, 
                                  fill = one.geno)) + 
  stat_dotsinterval(side = "bottom", 
                    slab_linewidth = NA) +
  scale_fill_manual(values = c("#003566",
                               "#ada7c9")) +
  labs(title = "Expression variation of genes with accessible and non-accessible promoters",
       x = "Accessible promoter (≤ 1 kb)") +
  theme_pubr() + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("TRUE", "FALSE")),
                     label = "p.format",
                     step.increase = 0.05) +
  theme(legend.position = "none") 



# Look at relationship between distance to TSS and baseMean/lfcSE
### Check first if the distanceToTSS varies for each gene between the genotypes
geno.1.promoter$geno <- as.factor("1") 
geno.2.promoter$geno <- as.factor("2")
geno.3.promoter$geno <- as.factor("3") 
geno.4.promoter$geno <- as.factor("4")



# TopGO of peaks