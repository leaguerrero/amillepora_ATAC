# WGCNA of Chapter 2 followed by igraph to create hub scores
library(WGCNA)
library(DESeq2)
library(tidyr)
library(tidyverse)
library(broom)
library(purrr)
options(stringsAsFactors = FALSE)

#Enable multithread
enableWGCNAThreads()

# Load the expression data
# File Set-up
my.dir <- "~/Lab Notebook/Chapter2/Data_analysis/TagSeq/Gene_count_files/"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.table("/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)



# Create sample table
sampleNames <- my.metadata$ID
my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files, 
                             condition = my.metadata)
# Convert variables to factors
factorVars <- c("condition.genotype", "condition.ext.batch", "condition.stress.txt", "condition.acclimation.temp")
my.sampleTable[, factorVars] <- lapply(my.sampleTable[, factorVars], factor)
colnames(my.sampleTable)

# Create DESeqDataSet Object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable,
  directory = my.dir,
  design = ~ condition.stress.group*condition.genotype)
# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

amill.counts.matrix <- ddsHTSeq@assays@data@listData$counts

# Transform matrix
amill.counts.mat <- varianceStabilizingTransformation(amill.counts.matrix , blind = TRUE)

# Transpose the expression data
counts.df <- as.data.frame(amill.counts.mat)
counts.df <- tibble::rownames_to_column(counts.df, "Genes")
datExpr.0 = as.data.frame(t(counts.df[, -c(1)]))
names(datExpr.0) = counts.df$Genes;
rownames(datExpr.0) = names(counts.df)[-c(1)]
sample.names <- names(counts.df)[-c(1)]

gsg.1 = goodSamplesGenes(datExpr.0, verbose = 3)
gsg.1$allOK
rm(gsg.1)
# Check Sample Outliters
sampleTree.1 = hclust(dist(datExpr.0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
plot(sampleTree.1, main = "Sample clustering to detect outliers",
     sub="",
     xlab="",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red")
# dev.off()
# # No outliers
datExpr <- datExpr.0
# 

# 
# Picking beta for expression data
# Choose a set of soft-thresholding powers
powers.e = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft.e = pickSoftThreshold(datExpr, powerVector = powers.e, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft.e$fitIndices[,1], -sign(sft.e$fitIndices[,3])*sft.e$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft.e$fitIndices[,1], -sign(sft.e$fitIndices[,3])*sft.e$fitIndices[,2],
     labels=powers.e,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.82,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft.e$fitIndices[,1], sft.e$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.e$fitIndices[,1], sft.e$fitIndices[,5], labels=powers.e, cex=cex1,col="red")
abline(h=0.90,col="red")

# # based on plot, I will pick a threshold of 2
# 
# Constructing the gene network
ge.network= blockwiseModules(datExpr, power = 2,
                             TOMType = "signed", minModuleSize = 100,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = FALSE,
                             verbose = 3)

# Restart R Session
# library(WGCNA)
# search(): ".GlobalEnv", "package:WGCNA"
table(ge.network$colors)


sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(ge.network$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ge.network$dendrograms[[1]], mergedColors[ge.network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels.e = ge.network$colors
moduleColors.e = labels2colors(ge.network$colors)
ME.es = ge.network$MEs;
geneTree.e = ge.network$dendrograms[[1]]


softPower = 2 #Chosen in the graphs before
adjacency = adjacency(datExpr, power = softPower, type = "signed")
#Transforming the adjacency matrix in a topological overlap
TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
 rm(adjacency)

# ID Heat response modules
 
 #load igraph
 library(igraph)
graph <- igraph::graph_from_adjacency_matrix(TOM) 
plot(graph)

hub.scores <- hub_score(graph, scale = TRUE, weights = NULL, options = arpack_defaults())
