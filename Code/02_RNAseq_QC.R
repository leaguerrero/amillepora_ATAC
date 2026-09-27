# 02_RNAseq_QC.R
#
# Purpose:
#   Perform exploratory quality-control analyses of the RNA-seq
#   dataset before differential-expression analysis.
#
#   This script:
#     1. Reads the filtered RNA-seq count matrix from Script 01
#     2. Calculates DESeq2-normalized counts
#     3. Performs variance-stabilizing transformation (VST)
#     4. Performs whole-transcriptome PCA
#     5. Partitions gene-expression variance among experimental
#        variables using variancePartition
#
# Inputs from Script 01:
#   - RNAseq_filtered_count_matrix.csv
#   - RNAseq_sample_table.csv
#
# Outputs:
#   - RNAseq_normalized_counts.csv
#   - RNAseq_VST_counts.csv
#   - RNAseq_PCA_heat_stress_scores.csv
#   - RNAseq_PCA_heat_stress.pdf
#   - RNAseq_PCA_heat_stress.tiff
#   - RNAseq_variance_partition_all_genes.csv
#   - RNAseq_variance_partition_mean_percent.csv
#   - RNAseq_variance_partition.pdf
#   - RNAseq_variance_partition.tiff
#
# Expected variance-partition results reported in manuscript:
#   Genotype:              ~22.3%
#   Acute heat stress:     ~17.6%
#   Heat-stress assay day:  ~9.1%
#   Extraction batch:       ~6.8%
#   Thermal history:        ~5.0%



# 1. Load packages

library(DESeq2)
library(variancePartition)
library(ggplot2)



# 2. Define input and output directories


project.dir <- "./Code_and_data/"

count.matrix.file <- file.path(
  project.dir,
  "./Data/01_count_matrix",
  "RNAseq_filtered_count_matrix.csv"
)

sample.table.file <- file.path(
  project.dir,
  "./Data/01_count_matrix",
  "RNAseq_sample_table.csv"
)

output.dir <- file.path(
  project.dir,
  "./Data/02_QC"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Read filtered count matrix


filtered.counts <- read.csv(
  count.matrix.file,
  check.names = FALSE
)

cat(
  "Dimensions of filtered count file:",
  nrow(filtered.counts),
  "genes x",
  ncol(filtered.counts) - 1,
  "samples\n"
)

# Dimensions of filtered count file: 19649 genes x 58 samples

# The first column contains gene IDs

gene.ids <- filtered.counts$gene


# Remove the gene-ID column to create a numeric count matrix

count.matrix <- as.matrix(
  filtered.counts[, -1]
)

rownames(count.matrix) <- gene.ids


# HTSeq counts should be whole numbers

storage.mode(count.matrix) <- "integer"



# 4. Read sample table

sample.table <- read.csv(
  sample.table.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Number of samples in sample table:",
  nrow(sample.table),
  "\n"
)

# Number of samples in sample table: 58 

cat(
  "\nSample-table column names:\n"
)

print(
  colnames(sample.table)
)

# [1] "sampleName"        "fileName"          "acclimation.temp" 
# [4] "acclimation.class" "acclimation.txt"   "stress.txt"       
# [7] "stress.group"      "stress.day"        "ext.date"         
# [10] "ext.batch"         "genotype" 


# 5. Confirm that count matrix and sample table agree


if (ncol(count.matrix) != nrow(sample.table)) {
  stop(
    "Number of count-matrix columns does not match number of sample-table rows."
  )
}


# Reorder the sample table to match the count matrix

sample.table <- sample.table[
  match(
    colnames(count.matrix),
    sample.table$sampleName
  ),
]


# Confirm that every count-matrix sample was found

if (any(is.na(sample.table$sampleName))) {
  stop(
    "At least one count-matrix sample was not found in the sample table."
  )
}


# Confirm identical sample order

if (!all(sample.table$sampleName == colnames(count.matrix))) {
  stop(
    "Sample order does not match between count matrix and sample table."
  )
}


rownames(sample.table) <- sample.table$sampleName



# 6. Create clearly named experimental variables

# These correspond to the variables used in the manuscript:
#
#   genotype
#   heat-stress assay day
#   extraction batch
#   thermal-history treatment
#   acute heat-stress exposure
#
# The original metadata column names are retained in the
# sample table. These new columns simply provide clearer names
# for the statistical analyses below.


sample.table$genotype <- factor(
  sample.table$genotype
)

sample.table$stress_day <- factor(
  sample.table$stress.day
)

sample.table$extraction_batch <- factor(
  sample.table$ext.batch
)

sample.table$thermal_history <- factor(
  sample.table$acclimation.txt
)

sample.table$heat_stress <- factor(
  sample.table$stress.txt
)


# Inspect the levels of each variable

cat(
  "\nGenotype levels:\n"
)
print(
  levels(sample.table$genotype)
)

# [1] "1" "2" "3" "4" "5" "6" "7" "8"

cat(
  "\nHeat-stress assay day levels:\n"
)
print(
  levels(sample.table$stress_day)
)

# [1] "D00" "D07" "D14" "D21"


cat(
  "\nExtraction-batch levels:\n"
)
print(
  levels(sample.table$extraction_batch)
)

# [1] "1" "2" "3" "4" "5"

cat(
  "\nThermal-history levels:\n"
)
print(
  levels(sample.table$thermal_history)
)

# [1] "7DA_0DB"  "7DA_14DB" "7DA_7DB"  "control" 

cat(
  "\nHeat-stress treatment levels:\n"
)
print(
  levels(sample.table$heat_stress)
)
# [1] "24.5" "36"  


# 7. Construct a DESeq2 object for normalization and QC

# No differential-expression model is fit here.
# The design is ~ 1.


dds.qc <- DESeqDataSetFromMatrix(
  countData = count.matrix,
  colData = sample.table,
  design = ~ 1
)

# 8. Estimate DESeq2 size factors


dds.qc <- estimateSizeFactors(
  dds.qc
)


# Inspect size factors

cat(
  "\nDESeq2 size factors:\n"
)

print(
  sizeFactors(dds.qc)
)

# 01_27_Y   01_28_R   01_40_Y   01_42_R   01_45_Y   01_46_Y   01_51_Y   01_52_Y 
# 1.2998584 0.3825635 0.9137377 0.7259135 2.3225146 0.9524503 2.4615931 0.5068095 
# 01_54_Y   01_58_Y   01_63_Y   01_68_Y   01_78_R   01_92_Y   02_29_Y   02_34_Y 
# 0.7828314 0.7193197 0.5810685 0.8115313 0.7454714 1.0489890 1.0441529 2.5724126 
# 02_35_R   02_47_Y   02_65_R   02_79_R   03_29_R   03_48_Y   03_62_R   03_72_R 
# 1.7460645 0.8895893 0.8949021 1.1564632 2.6652877 0.4791176 2.3447350 0.4065928 
# 03_90_Y   03_91_R   04_48_R   04_65_Y   04_70_Y   04_95_Y   05_43_R   05_46_R 
# 1.0557439 1.3437273 0.8589865 0.7703772 1.6341477 1.7935754 0.6754654 2.0488862 
# 05_50_Y   05_62_Y   05_82_Y   05_88_R   06_22_Y   06_32_R   06_59_R   06_64_R 
# 0.4207898 0.9817953 2.0822704 0.3739338 1.4141079 0.9946759 0.7877487 0.3157539 
# 06_67_R   06_72_Y   06_82_R   06_92_R   07_31_R   07_36_R   07_39_R   07_53_Y 
# 1.0973007 1.3095512 2.3540903 2.8025500 3.0118909 1.0432914 0.8515427 0.4438296 
# 07_64_Y   07_88_Y   08_19_R   08_30_R   08_33_Y   08_36_Y   08_45_R   08_50_R 
# 0.4199774 0.9686339 0.7321219 2.0741506 1.9684520 0.6494242 0.9793374 1.2209599 
# 08_85_Y   08_95_R 
# 0.8675558 1.5426242 

# 9. Extract normalized counts

normalized.counts <- counts(
  dds.qc,
  normalized = TRUE
)


cat(
  "\nNormalized-count matrix dimensions:",
  nrow(normalized.counts),
  "genes x",
  ncol(normalized.counts),
  "samples\n"
)

# Normalized-count matrix dimensions: 19649 genes x 58 samples

# Save normalized counts

normalized.counts.output <- data.frame(
  gene = rownames(normalized.counts),
  normalized.counts,
  check.names = FALSE
)

write.csv(
  normalized.counts.output,
  file = file.path(
    output.dir,
    "RNAseq_normalized_counts.csv"
  ),
  row.names = FALSE
)

# PART A: Whole-transcriptome PCA
# 10. Variance-stabilizing transformation

# VST is used for PCA so that genes with very large counts do
# not dominate the analysis solely because of count magnitude.


vst.data <- vst(
  dds.qc,
  blind = TRUE
)

vst.matrix <- assay(
  vst.data
)


cat(
  "\nVST matrix dimensions:",
  nrow(vst.matrix),
  "genes x",
  ncol(vst.matrix),
  "samples\n"
)

# VST matrix dimensions: 19649 genes x 58 samples

# Save VST-transformed matrix

vst.output <- data.frame(
  gene = rownames(vst.matrix),
  vst.matrix,
  check.names = FALSE
)

write.csv(
  vst.output,
  file = file.path(
    output.dir,
    "RNAseq_VST_counts.csv"
  ),
  row.names = FALSE
)


# 11. Perform PCA

# Samples are colored by acute heat-stress treatment.

pca.data <- plotPCA(
  vst.data,
  intgroup = "heat_stress",
  returnData = TRUE
)


# Extract percent variance explained by PC1 and PC2

percent.variance <- round(
  100 * attr(
    pca.data,
    "percentVar"
  ),
  1
)


cat(
  "\nVariance explained by PC1:",
  percent.variance[1],
  "%\n"
)

# Variance explained by PC1: 50.2 %

cat(
  "Variance explained by PC2:",
  percent.variance[2],
  "%\n"
)

# Variance explained by PC2: 12.1 %

# Save PCA coordinates

write.csv(
  pca.data,
  file = file.path(
    output.dir,
    "RNAseq_PCA_heat_stress_scores.csv"
  ),
  row.names = TRUE
)

# 12. Plot PCA


pca.plot <- ggplot(
  pca.data,
  aes(
    x = PC1,
    y = PC2,
    color = heat_stress
  )
) +
  geom_point(
    size = 4,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "24.5" = "#bfd7ea",
      "36" = "#ff5a5f"
    ),
    labels = c(
      "24.5" = "Control",
      "36" = "Heat stress"
    )
  ) +
  labs(
    title = "PCA of A. millepora gene expression",
    x = paste0(
      "PC1 (",
      percent.variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent.variance[2],
      "%)"
    ),
    color = "Heat-stress treatment"
  ) +
  theme_classic(
    base_size = 16
  )

pca.plot


# Save PCA as PDF

ggsave(
  filename = file.path(
    output.dir,
    "RNAseq_PCA_heat_stress.pdf"
  ),
  plot = pca.plot,
  width = 8,
  height = 6,
  units = "in"
)


# Save PCA as TIFF

ggsave(
  filename = file.path(
    output.dir,
    "RNAseq_PCA_heat_stress.tiff"
  ),
  plot = pca.plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)


# PART B: Variance partitioning
# 13. Log2-transform DESeq2-normalized counts

# This reproduces the expression values used in the original
# variancePartition analysis.


log.normalized.counts <- log2(
  normalized.counts + 1
)


# 14. Define variance-partition model


variance.formula <- ~
  genotype +
  stress_day +
  extraction_batch +
  thermal_history +
  heat_stress


cat(
  "\nVariance-partition formula:\n"
)

print(
  variance.formula
)

# ~genotype + stress_day + extraction_batch + thermal_history + 
  # heat_stress

# 15. Run variancePartition

# This estimates the fraction of expression variance associated
# with each experimental variable separately for each gene.


variance.results <- fitExtractVarPartModel(
  exprObj = log.normalized.counts,
  formula = variance.formula,
  data = sample.table
)


# 16. Calculate mean variance explained across genes


mean.variance <- colMeans(
  variance.results
) * 100


mean.variance.table <- data.frame(
  variable = names(mean.variance),
  mean_percent_variance = as.numeric(mean.variance),
  row.names = NULL
)


cat(
  "\nMean percent variance explained across genes:\n"
)

print(
  mean.variance.table
)
# variable mean_percent_variance
# 1         genotype             22.273018
# 2       stress_day              9.140355
# 3 extraction_batch              6.830386
# 4  thermal_history              4.966133
# 5      heat_stress             17.635454
# 6        Residuals             39.154654


# 17. Save variance-partition results


variance.results.output <- data.frame(
  gene = rownames(variance.results),
  variance.results,
  check.names = FALSE
)

write.csv(
  variance.results.output,
  file = file.path(
    output.dir,
    "RNAseq_variance_partition_all_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  mean.variance.table,
  file = file.path(
    output.dir,
    "RNAseq_variance_partition_mean_percent.csv"
  ),
  row.names = FALSE
)



# 18. Plot variance-partition results


variance.plot <- plotVarPart(
  variance.results
) +
  labs(
    title = "Variance partition of gene expression",
    x = NULL,
    y = "Fraction of variance explained"
  ) +
  theme_bw(
    base_size = 14
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.minor = element_blank()
  )

variance.plot


# Save variance-partition plot as PDF

ggsave(
  filename = file.path(
    output.dir,
    "RNAseq_variance_partition.pdf"
  ),
  plot = variance.plot,
  width = 8,
  height = 6,
  units = "in"
)


# Save variance-partition plot as TIFF

ggsave(
  filename = file.path(
    output.dir,
    "RNAseq_variance_partition.tiff"
  ),
  plot = variance.plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)



# 19. Final checks


cat(
  "RNA-seq QC complete\n"
)


cat(
  "Samples:",
  ncol(count.matrix),
  "\n"
)

cat(
  "Genes:",
  nrow(count.matrix),
  "\n\n"
)

cat(
  "Mean variance explained:\n"
)

print(
  mean.variance.table
)



# 20. Record R and package versions

capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "RNAseq_QC_sessionInfo.txt"
  )
)
