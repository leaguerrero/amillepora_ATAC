# ATAC Seq and Gene expression plasticity TOPGO

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
library(patchwork)
library(introdataviz)
library(topGO)
library(tidyverse)

# Load in data sets
# Transcripts
gff.file <- "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amil.all.maker.noseq.gff"
txdb <- makeTxDbFromGFF(gff.file, format = "gff3") 
a.millepora.txs <- transcripts(txdb)

# Consensus peaks
library(DiffBind)
samples <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/ATAC/sampleSheet.csv", header = TRUE)
samples
DBdata <- dba(sampleSheet=samples)
DBdata

DBdata <- dba.count(DBdata,  bUseSummarizeOverlaps=TRUE) # normalized counts per https://www.biostars.org/p/401188/
consensus <- dba(DBdata, DBdata$masks$Consensus)
consensus.peaks <- dba.peakset(DBdata, minOverlap=2, peak.caller = "macs",peak.format = "macs",  bRetrieve=TRUE) # produces GR object of consensus peaks
consensus.peaks.anno <- annotatePeak(consensus.peaks, TxDb = txdb)
consensus.peaks.anno.df <- as.data.frame(consensus.peaks.anno@anno@elementMetadata@listData)
consensus.peaks.anno.df$gene <- consensus.peaks.anno.df$geneId

# Consensus promoters
consensus.promoters <- consensus.peaks.anno.df[(grepl("promoter", consensus.peaks.anno.df$annotation, ignore.case = TRUE)),] #3973


# Gene Expression Preparation
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

# Extract normalized counts matrix
# Estimate size factors
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)

# Estimate normalization factors
ddsHTSeq <- estimateDispersions(ddsHTSeq)
normalized.counts <- as.data.frame(counts(ddsHTSeq, normalized = TRUE))
normalized.counts$gene <- rownames(normalized.counts) 
normalized.counts$standard.deviation <- apply(normalized.counts[,1:58], 1, sd)

# ID HRGs
dds <- DESeq(ddsHTSeq)

# Save gene set for topGO analysis
total.set.gene.names <- rownames(assay(dds))
total.set.gene.names.df <- data.frame(gene = total.set.gene.names)
# write.csv(total.set.gene.names.df, file = '~/Lab Notebook/Chapter2/Data_analysis/Ch2.total.set.gene.names.csv', row.names = FALSE)


# Select for heat response genes
hs.results <- lfcShrink(dds,
                        contrast = c("condition.stress.txt", "36", "24.5"),
                        type = "ashr")

heat.response.results <- hs.results[complete.cases(hs.results),]
heat.response.results$gene <- rownames(heat.response.results)

in.peaks.two.genos.base.mean<- heat.response.results$gene %in% geneID_2_outta_4 #  genes overlap
heat.response.results$two.genos <- in.peaks.two.genos.base.mean
sum(heat.response.results$two.genos) # 896

# Filter out rows with high adj p-value
hist(heat.response.results$padj)
padj.filtered.heat.response.results <- heat.response.results[heat.response.results$padj < 0.05, ]
hist(padj.filtered.heat.response.results$padj)

# Filter for change in heat 
amillepora.HRGs <- subset(padj.filtered.heat.response.results,
                      abs(log2FoldChange) >= 2)

# Uncertainty over background gene set to use: DC with Rachael

# Gene Ontology (GO) terms for heat response genes

# Set-up
# Read in annotation table and reformat GO column
annos <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",sep="\t",IDsep=",")

# Read in set of all genes and make 'background' gene set
all <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/Ch2.total.set.gene.names.csv")
allgenes <- all$gene

# Modify heat response gene set for topGO analysis

amillepora.hrgIG <- factor(as.numeric(allgenes%in%amillepora.HRGs$gene))
names(amillepora.hrgIG) <- allgenes

##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.hrgIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hr.genes.GO_BP.csv",sep=""))

##MF
GOamillepora <- new("topGOdata",ontology="MF",allGenes=amillepora.hrgIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hr.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.hrgIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hr.genes.GO_CC.csv",sep=""))


## TopGO of peaks (2 out of 4 genos)

# Create background set of genes with peaks with DiffBind (uses DESeq2)
## update the backgroung gene set to be all genes:
library(DiffBind)
library(ChIPpeakAnno)

samples <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/ATAC/sampleSheet.csv")
samples

# DBdata <- dba(sampleSheet=samples)
# plot(DBdata)
# DBdata <- dba.count(DBdata)
# consensus <- dba(DBdata, DBdata$masks$Consensus)
# consensus.peaks <- dba.peakset(consensus, bRetrieve=TRUE) # produces GR object of consensus peaks

# Load in transcripts
gff.file <- "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amil.all.maker.noseq.gff"
txdb <- makeTxDbFromGFF(gff.file, format = "gff3") 
a.millepora.txs <- transcripts(txdb)

# consensus.peaks.anno <- annotatePeak(consensus.peaks, TxDb = txdb) # background set of genes in peaks
peaks.allgenes <- a.millepora.txs$tx_name
peaks.allgenes <- sub("-RA$", "", peaks.allgenes)
# peaks in promoters in at least 2 out 4 individuals
head(geneID_2_outta_4.df)
amillepora.promoter.peaksIG <- factor(as.numeric(peaks.allgenes%in%geneID_2_outta_4.df$gene))
names(amillepora.promoter.peaksIG) <- peaks.allgenes

##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
#write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.promoter.genes.GO_BP.csv",sep=""))

##MF
GOamillepora <- new("topGOdata",ontology="MF",amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.promoter.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.promoter.genes.GO_CC.csv",sep=""))

## TopGo of HRGs with promoter peaks
annos <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",sep="\t",IDsep=",")

# Read in set of all genes and make 'background' gene set
all <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/Ch2.total.set.gene.names.csv")
allgenes <- all$gene

# Gene set that has a large response to heat stress and have accessible promoters
hs.results <- lfcShrink(dds,
                        contrast = c("condition.stress.txt", "36", "24.5"),
                        type = "ashr")

heat.response.results <- hs.results[complete.cases(hs.results),]
heat.response.results$gene <- rownames(heat.response.results)

# Merge this with all genes in peaks (in at least 2 out of 4 samples)
# ID genes in peaks in genes 

### Insert other code to get to this point here ###
in.peaks.two.genos.base.mean<- heat.response.results$gene %in% geneID_2_outta_4 #  genes overlap
heat.response.results$two.genos <- in.peaks.two.genos.base.mean
sum(heat.response.results$two.genos) # 896
padj.filtered.heat.response.results <- heat.response.results[heat.response.results$padj < 0.05, ]
padj.filtered.heat.response.results
sum(padj.filtered.heat.response.results$two.genos) #429
genes.hrg.peaks <- padj.filtered.heat.response.results$gene[padj.filtered.heat.response.results$two.genos == TRUE]

amillepora.hrg.peaks.IG <- factor(as.numeric(allgenes%in%genes.hrg.peaks))
sum(amillepora.hrg.peaks.IG == 1)
names(amillepora.hrg.peaks.IG ) <- allgenes

##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.hrg.peaks.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hrg.peaks.promoter.genes.GO_BP.csv",sep=""))

##MF
GOamillepora <- new("topGOdata",ontology="MF",amillepora.hrg.peaks.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hrg.peaks.promoter.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.hrg.peaks.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hrg.peaks.promoter.genes.GO_CC.csv",sep=""))


## What about all genic peaks?
# First, isolate peaks id'ed as genic
geno.1.genic.peak <- filter(peaks.df.1, genic == "TRUE") 
geno.2.genic.peak <- filter(peaks.df.2, genic == "TRUE")  
geno.3.genic.peak <- filter(peaks.df.3, genic == "TRUE") 
geno.4.genic.peak <- filter(peaks.df.4, genic == "TRUE")

# Filter genic peaks for unique gene names and adjusted pValue

pval.threshold <- -log10(0.00005)

geno.1.genic.peak <- filter(geno.1.genic.peak, pValue >= pval.threshold) 
genes.geno.1 <- unique(data.frame(gene = geno.1.genic.peak$geneId))
genes.geno.1$geno<- as.factor("1")
dim(genes.geno.1)

geno.2.genic.peak<- filter(geno.2.genic.peak, pValue >= pval.threshold)
genes.geno.2 <- unique(data.frame(gene = geno.2.genic.peak$geneId))
genes.geno.2$geno<- as.factor("2")
dim(genes.geno.2)

geno.3.genic.peak <- filter(geno.3.genic.peak, pValue >= pval.threshold) 
genes.geno.3 <- unique(data.frame(gene = geno.3.genic.peak$geneId))
genes.geno.3$geno<- as.factor("3")
dim(genes.geno.3)

geno.4.genic.peak <- filter(geno.4.genic.peak, pValue >= pval.threshold) 
genes.geno.4 <- unique(data.frame(gene = geno.4.genic.peak$geneId))
genes.geno.4$geno<- as.factor("4")
dim(genes.geno.4) 

# ID Peaks in 2 out of 4 genotypes
geneID_geno1 <- unique(as.data.frame(geno.1.genic.peak$geneId))
colnames(geneID_geno1) <- "geneID"
geneID_geno2 <- unique(as.data.frame(geno.2.genic.peak$geneId))
colnames(geneID_geno2) <- "geneID"
geneID_geno3 <- unique(as.data.frame(geno.3.genic.peak$geneId))
colnames(geneID_geno3) <- "geneID"
geneID_geno4 <- unique(as.data.frame(geno.4.genic.peak$geneId))
colnames(geneID_geno4) <- "geneID"

geneID_1_2 <- dplyr::inner_join(geneID_geno1, geneID_geno2, by = "geneID")
geneID_1_3 <- dplyr::inner_join(geneID_geno1, geneID_geno3, by = "geneID")
geneID_1_4 <- dplyr::inner_join(geneID_geno1, geneID_geno4, by = "geneID")

geneID_2_3 <- dplyr::inner_join(geneID_geno2, geneID_geno3, by = "geneID")
geneID_2_4 <- dplyr::inner_join(geneID_geno2, geneID_geno4, by = "geneID")

geneID_3_4 <- dplyr::inner_join(geneID_geno3, geneID_geno4, by = "geneID")

geneID_2_outta_4 <- rbind(geneID_1_2, 
                          geneID_1_3, 
                          geneID_1_4, 
                          geneID_2_3, 
                          geneID_2_4,
                          geneID_3_4)

geneID_2_outta_4 <- unique(geneID_2_outta_4$geneID) #2868

geneID_2_outta_4.df <- as.data.frame(geneID_2_outta_4)
colnames(geneID_2_outta_4.df) <- "gene"

# topGo of all genic peaks
# peaks in promoters in at least 2 out 4 individuals
head(geneID_2_outta_4.df)
amillepora.genic.peaksIG <- factor(as.numeric(peaks.allgenes%in%geneID_2_outta_4.df$gene))
names(amillepora.genic.peaksIG) <- peaks.allgenes

##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.genic.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.genes.GO_BP.csv",sep=""))

##MF
GOamillepora <- new("topGOdata",ontology="MF",amillepora.genic.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.genic.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.genes.GO_CC.csv",sep=""))

## TopGo of HRGs with genic peaks
annos <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")

# Combine genic peaks with gene expression genic results
in.peaks.two.genos.base.mean<- heat.response.results$gene %in% geneID_2_outta_4 #  genes overlap
heat.response.results$two.genos <- in.peaks.two.genos.base.mean
sum(heat.response.results$two.genos) # 2473

# Filter out rows with high adj p-value
hist(heat.response.results$padj)
padj.filtered.heat.response.results <- heat.response.results[heat.response.results$padj < 0.05, ]
hist(padj.filtered.heat.response.results$padj)
sum(padj.filtered.heat.response.results$two.genos) #1254

# TopGo
genes.hrg.peaks <- padj.filtered.heat.response.results$gene[padj.filtered.heat.response.results$two.genos == TRUE]

amillepora.hrg.genic.peaks.IG <- factor(as.numeric(allgenes%in%genes.hrg.peaks))
sum(amillepora.hrg.genic.peaks.IG == 1)
names(amillepora.hrg.peaks.IG ) <- allgenes


##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.genic.peaks.hrg.IG , nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.hr.genes.GO_BP.csv",sep=""))


##MF
GOamillepora <- new("topGOdata",ontology="MF",amillepora.genic.peaks.hrg.IG , nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
#write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.hr.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.genic.peaks.hrg.IG , nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
# write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.peaks.genic.hr.genes.GO_CC.csv",sep=""))

### Consensus peak promoters GO Terms
# Remove duplicate gene names
unique.consensus.promoters <- consensus.promoters[!duplicated(consensus.promoters$gene),]
unique.consensus.promoters.genes <- unique.consensus.promoters$gene

peaks.allgenes <- a.millepora.txs$tx_name
peaks.allgenes <- sub("-RA$", "", peaks.allgenes)
# unique consensus promoter peaks 
head(unique.consensus.promoters.genes)
amillepora.promoter.peaksIG <- factor(as.numeric(peaks.allgenes%in%unique.consensus.promoters$gene))
names(amillepora.promoter.peaksIG) <- peaks.allgenes

##BP
GOamillepora <- new("topGOdata",ontology="BP",allGenes=amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/Consensus peaks/Promoters/amillepora.peaks.promoter.genes.GO_BP.csv",sep=""))

##MF
GOamillepora <- new("topGOdata",ontology="MF",amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score), numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/Consensus peaks/Promoters/amillepora.peaks.promoter.genes.GO_MF.csv",sep=""))

##CC
GOamillepora <- new("topGOdata",ontology="CC",allGenes=amillepora.promoter.peaksIG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
amilleporaFisher <- getSigGroups(GOamillepora,test.stat)
amillepora.res <- GenTable(GOamillepora,classic=amilleporaFisher,topNodes=length(amilleporaFisher@score),numChar=100)
amillepora.filt <- amillepora.res[amillepora.res$classic<0.01 & amillepora.res$Significant>=10,]
write.csv(amillepora.filt,paste("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/Consensus peaks/Promoters/amillepora.peaks.promoter.genes.GO_CC.csv",sep=""))






