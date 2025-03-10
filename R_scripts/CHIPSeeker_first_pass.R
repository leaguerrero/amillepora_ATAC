## ATACseq
# CHIP_seeker
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)

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

##### Assessing each genotype individually
# ChIP peaks coverage plot
library(patchwork)
covplot1 <- covplot(peak.gr.list[[1]], 
        weightCol="V5", 
        chrs=c("chr1","chr2","chr3", "chr4","chr5", "chr6", "chr7", "chr8","chr9", "chr10", "chr11", "chr12","chr13", "chr14"),
        title = "Geno 1: ChIP Peaks over Chromosomes",
        fill_color = "#f72585")

covplot2 <- covplot(peak.gr.list[[2]], 
        weightCol="V5", 
        chrs=c("chr1","chr2","chr3", "chr4","chr5", "chr6", "chr7", "chr8","chr9", "chr10", "chr11", "chr12","chr13", "chr14"),
        title = "Geno 2: ChIP Peaks over Chromosomes",
        fill_color = "#7209b7")

covplot3 <- covplot(peak.gr.list[[3]], 
        weightCol="V5", 
        chrs=c("chr1","chr2","chr3", "chr4","chr5", "chr6", "chr7", "chr8","chr9", "chr10", "chr11", "chr12","chr13", "chr14"),
        title = "Geno 3: ChIP Peaks over Chromosomes",
        fill_color = "#3a0ca3")

covplot4 <- covplot(peak.gr.list[[4]], 
        weightCol="V5", 
        chrs=c("chr1","chr2","chr3", "chr4","chr5", "chr6", "chr7", "chr8","chr9", "chr10", "chr11", "chr12","chr13", "chr14"),
        title = "Geno 4: ChIP Peaks over Chromosomes",
        fill_color = "#4361ee")
(covplot1 |covplot2) /
  (covplot3|covplot4)

# Peak Annotation
peakAnno.1 <- annotatePeak(peak.gr.list[[1]], tssRegion=c(-3000, 3000),
                         TxDb=txdb)

peakAnno.2 <- annotatePeak(peak.gr.list[[2]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

peakAnno.3 <- annotatePeak(peak.gr.list[[3]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)

peakAnno.4 <- annotatePeak(peak.gr.list[[4]], tssRegion=c(-3000, 3000),
                           TxDb=txdb)


anno.plot.1 <- upsetplot(peakAnno.1, vennpie=TRUE)
anno.plot.2 <- upsetplot(peakAnno.2, vennpie=TRUE)
anno.plot.3 <- upsetplot(peakAnno.3, vennpie=TRUE)
anno.plot.4 <- upsetplot(peakAnno.4, vennpie=TRUE)

(anno.plot.1 |anno.plot.2) /
  (anno.plot.3|anno.plot.4)


((anno.plot.1 |anno.plot.2) + theme(plot.margin = unit(c(10,40,10,40), "pt"))) /
  ( (anno.plot.3|anno.plot.4) + theme(plot.margin = unit(c(10,40,10,40), "pt")))


## classify this as accessible across the whole gene, or by the TSS region
### Since with ATAC-Seq enrichment can be preferential to the TSS, 
### it's possible you want to rely on TSS enrichment as opposed to requiring 
### the entire gene to show accessibility before calling it open. 

peaks.df <- as.data.frame(peakAnno.1@anno@elementMetadata@listData)
genomic.annotation.detail.df <- peakAnno.1@detailGenomicAnnotation
peaks.df <- cbind(peaks.df, genomic.annotation.detail.df)

