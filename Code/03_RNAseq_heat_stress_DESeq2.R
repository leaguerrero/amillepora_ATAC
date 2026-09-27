# 03_RNAseq_heat_stress_DESeq2.R
#
# Purpose:
#   Identify genes that respond to acute heat stress in
#   Acropora millepora using the final DESeq2 model.
#
#   The model includes:
#     - genotype
#     - heat-stress assay day
#     - extraction batch
#     - thermal-history treatment
#     - acute heat-stress treatment
#
#   Acute heat stress is the primary effect of interest.
#
# Inputs from Script 01:
#   - RNAseq_filtered_count_matrix.csv
#   - RNAseq_sample_table.csv
#
# Outputs:
#   - RNAseq_DESeq2_sample_table.csv
#   - RNAseq_DESeq2_design_matrix.csv
#   - RNAseq_heat_stress_all_results.csv
#   - RNAseq_heat_stress_DEGs.csv
#   - RNAseq_heat_stress_model.rds
#   - RNAseq_heat_stress_model_diagnostics.pdf
#   - RNAseq_heat_stress_model_summary.txt
#
# Expected checkpoint:
#   - 19,649 genes tested
#   - 2,611 final heat-responsive genes
#
# Final DEG criteria:
#   adjusted p-value < 0.01
#   absolute log2 fold change > 2
#   dispersion outliers excluded

# 1. Load package
library(DESeq2)

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
  "./Data/03_DESeq2"
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


# Save gene IDs separately

gene.ids <- filtered.counts$gene


# Remove gene column to create numeric count matrix

count.matrix <- as.matrix(
  filtered.counts[, -1]
)

rownames(count.matrix) <- gene.ids


# DESeq2 requires integer count data

storage.mode(count.matrix) <- "integer"


cat(
  "Count matrix:",
  nrow(count.matrix),
  "genes x",
  ncol(count.matrix),
  "samples\n"
)

# Count matrix: 19649 genes x 58 samples


# 4. Read sample metadata


sample.table <- read.csv(
  sample.table.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Sample table:",
  nrow(sample.table),
  "samples\n"
)

# Sample table: 58 samples


# 5. Match sample-table order to count-matrix order

sample.table <- sample.table[
  match(
    colnames(count.matrix),
    sample.table$sampleName
  ),
]


# Confirm every sample was found

if (any(is.na(sample.table$sampleName))) {
  stop(
    "At least one count-matrix sample was not found in the sample table."
  )
}


# Confirm sample order is identical

if (!all(sample.table$sampleName == colnames(count.matrix))) {
  stop(
    "Sample order does not match between count matrix and sample table."
  )
}


rownames(sample.table) <- sample.table$sampleName


cat(
  "Count matrix and sample table are in the same order.\n"
)


# 6. Create clearly named model variables


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



# 7. Set reference levels
# Thermal-history control is the reference for thermal history.
#
# 24.5 degrees C is the reference for acute heat stress,
# so the heat-stress coefficient represents:
#
#     36 degrees C versus 24.5 degrees C


sample.table$thermal_history <- relevel(
  sample.table$thermal_history,
  ref = "control"
)

sample.table$heat_stress <- relevel(
  sample.table$heat_stress,
  ref = "24.5"
)


# Inspect factor levels

cat(
  "\nThermal-history levels:\n"
)

print(
  levels(sample.table$thermal_history)
)
# [1] "control"  "7DA_0DB"  "7DA_14DB" "7DA_7DB" 

cat(
  "\nHeat-stress levels:\n"
)

print(
  levels(sample.table$heat_stress)
)
# [1] "24.5" "36" 


# 8. Inspect experimental structure


cat(
  "\nThermal history x heat stress:\n"
)

print(
  table(
    sample.table$thermal_history,
    sample.table$heat_stress
  )
)

# 24.5 36
# control    17 17
# 7DA_0DB     4  4
# 7DA_14DB    4  4
# 7DA_7DB     4  4

cat(
  "\nAssay day x heat stress:\n"
)

print(
  table(
    sample.table$stress_day,
    sample.table$heat_stress
  )
)

# 24.5 36
# D00    5  5
# D07    8  8
# D14    8  8
# D21    8  8

cat(
  "\nGenotype x heat stress:\n"
)

print(
  table(
    sample.table$genotype,
    sample.table$heat_stress
  )
)

# 24.5 36
# 1    7  7
# 2    3  3
# 3    3  3
# 4    2  2
# 5    3  3
# 6    4  4
# 7    3  3
# 8    4  4


# 9. Define final DESeq2 model


deseq.formula <- ~
  genotype +
  stress_day +
  extraction_batch +
  thermal_history +
  heat_stress


cat(
  "\nFinal DESeq2 model:\n"
)

print(
  deseq.formula
)
# ~genotype + stress_day + extraction_batch + thermal_history + 
#   heat_stress


# 10. Confirm that the design matrix is full rank


design.matrix <- model.matrix(
  deseq.formula,
  data = sample.table
)


cat(
  "\nDesign-matrix rank:",
  qr(design.matrix)$rank,
  "\n"
)

cat(
  "Number of design-matrix columns:",
  ncol(design.matrix),
  "\n"
)


if (qr(design.matrix)$rank != ncol(design.matrix)) {
  stop(
    "The DESeq2 design matrix is not full rank."
  )
}


cat(
  "Design matrix is full rank.\n"
)
# Design matrix is full rank.

# Save design matrix for documentation

design.matrix.output <- data.frame(
  sample = rownames(design.matrix),
  design.matrix,
  check.names = FALSE
)

write.csv(
  design.matrix.output,
  file = file.path(
    output.dir,
    "RNAseq_DESeq2_design_matrix.csv"
  ),
  row.names = FALSE
)



# 11. Confirm that the count filtering from Script 01 is intact

# Every gene should have at least 10 total reads across all
# samples because this filtering was performed in Script 01.

gene.total.counts <- rowSums(
  count.matrix
)


cat(
  "\nMinimum total count among retained genes:",
  min(gene.total.counts),
  "\n"
)

# Minimum total count among retained genes: 10 

if (any(gene.total.counts < 10)) {
  stop(
    "Some genes have fewer than 10 total counts. Check Script 01 output."
  )
}


cat(
  "All genes pass the >=10 total-count filter.\n"
)

# All genes pass the >=10 total-count filter.

# 12. Construct DESeq2 dataset


dds.heat <- DESeqDataSetFromMatrix(
  countData = count.matrix,
  colData = sample.table,
  design = deseq.formula
)


cat(
  "\nGenes included in DESeq2 analysis:",
  nrow(dds.heat),
  "\n"
)

# Genes included in DESeq2 analysis: 19649 

cat(
  "Samples included in DESeq2 analysis:",
  ncol(dds.heat),
  "\n"
)
# Samples included in DESeq2 analysis: 58 


# 13. Run DESeq2


dds.heat <- DESeq(
  dds.heat
)



# 14. Check coefficient convergence

# The original analysis had genes that did not initially
# converge. These were re-fit with a larger maximum number
# of Wald-test iterations.


cat(
  "\nInitial coefficient convergence:\n"
)

print(
  table(
    mcols(dds.heat)$betaConv,
    useNA = "ifany"
  )
)

# FALSE  TRUE 
# 294 19355 

failed.convergence <- sum(
  mcols(dds.heat)$betaConv == FALSE,
  na.rm = TRUE
)


cat(
  "Genes that did not initially converge:",
  failed.convergence,
  "\n"
)



# 15. Re-run Wald test if any genes failed to converge


if (failed.convergence > 0) {
  
  cat(
    "\nRe-running Wald test with maxit = 1000.\n"
  )
  
  dds.heat <- nbinomWaldTest(
    dds.heat,
    maxit = 1000
  )
}


# Check convergence again

cat(
  "\nFinal coefficient convergence:\n"
)

print(
  table(
    mcols(dds.heat)$betaConv,
    useNA = "ifany"
  )
)


final.failed.convergence <- sum(
  mcols(dds.heat)$betaConv == FALSE,
  na.rm = TRUE
)

# FALSE  TRUE 
# 154 19495 

cat(
  "Genes that failed to converge after final model fit:",
  final.failed.convergence,
  "\n"
)

# Genes that failed to converge after final model fit: 154 

# 16. Inspect DESeq2 coefficient names


cat(
  "\nDESeq2 coefficient names:\n"
)

print(
  resultsNames(dds.heat)
)
# [1] "Intercept"                           "genotype_2_vs_1"                    
# [3] "genotype_3_vs_1"                     "genotype_4_vs_1"                    
# [5] "genotype_5_vs_1"                     "genotype_6_vs_1"                    
# [7] "genotype_7_vs_1"                     "genotype_8_vs_1"                    
# [9] "stress_day_D07_vs_D00"               "stress_day_D14_vs_D00"              
# [11] "stress_day_D21_vs_D00"               "extraction_batch_2_vs_1"            
# [13] "extraction_batch_3_vs_1"             "extraction_batch_4_vs_1"            
# [15] "extraction_batch_5_vs_1"             "thermal_history_7DA_0DB_vs_control" 
# [17] "thermal_history_7DA_14DB_vs_control" "thermal_history_7DA_7DB_vs_control" 
# [19] "heat_stress_36_vs_24.5"             


# 17. Extract acute heat-stress results

# Contrast:
#     36 degrees C versus 24.5 degrees C
#
# Positive log2 fold change:
#     higher expression under heat stress
#
# Negative log2 fold change:
#     lower expression under heat stress
#
# The original analysis extracted DESeq2 results using
# alpha = 0.05. The final significance criterion is applied
# explicitly below as adjusted p-value < 0.01.

res.heat <- results(
  dds.heat,
  contrast = c(
    "heat_stress",
    "36",
    "24.5"
  ),
  alpha = 0.05
)


cat(
  "\nDESeq2 heat-stress result summary:\n"
)

print(
  summary(res.heat)
)

# out of 19649 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3564, 18%
# LFC < 0 (down)     : 5936, 30%
# outliers [1]       : 0, 0%
# low counts [2]     : 3810, 19%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# 18. Identify dispersion outliers

dispersion.outlier <- mcols(dds.heat)$dispOutlier %in% TRUE


high.disp.genes <- rownames(dds.heat)[
  dispersion.outlier
]


cat(
  "\nNumber of dispersion outlier genes:",
  length(high.disp.genes),
  "\n"
)

# Number of dispersion outlier genes: 362 


# 19. Create complete DESeq2 results table

res.heat.table <- as.data.frame(
  res.heat
)


res.heat.table$gene <- rownames(
  res.heat.table
)


res.heat.table$dispersionOutlier <- (
  res.heat.table$gene %in% high.disp.genes
)


# Put gene ID first

res.heat.table <- res.heat.table[
  ,
  c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "dispersionOutlier"
  )
]


# Sort by adjusted p-value

res.heat.table <- res.heat.table[
  order(
    res.heat.table$padj,
    na.last = TRUE
  ),
]



# 20. Identify significant heat-responsive genes
# Final criteria reported in the manuscript:
#
#   adjusted p-value < 0.01
#   absolute log2 fold change > 2
#   dispersion outliers excluded

heat.sig <- res.heat.table[
  !is.na(res.heat.table$padj) &
    res.heat.table$padj < 0.01 &
    abs(res.heat.table$log2FoldChange) > 2 &
    res.heat.table$dispersionOutlier == FALSE,
]


# Sort significant genes by adjusted p-value

heat.sig <- heat.sig[
  order(heat.sig$padj),
]

# 21. Summarize final heat-responsive gene set


number.heat.genes <- nrow(
  heat.sig
)


number.upregulated <- sum(
  heat.sig$log2FoldChange > 2
)


number.downregulated <- sum(
  heat.sig$log2FoldChange < -2
)


cat(
  "Final heat-responsive gene set\n"
)


cat(
  "Total heat-responsive genes:",
  number.heat.genes,
  "\n"
)

# Total heat-responsive genes: 2611 

cat(
  "Upregulated:",
  number.upregulated,
  "\n"
)

# Upregulated: 1170 

cat(
  "Downregulated:",
  number.downregulated,
  "\n"
)

# Downregulated: 1441 


# 22. Check whether dispersion outliers would otherwise
#     have passed the significance thresholds


significant.before.outlier.removal <- res.heat.table[
  !is.na(res.heat.table$padj) &
    res.heat.table$padj < 0.01 &
    abs(res.heat.table$log2FoldChange) > 2,
]


outliers.among.significant <- significant.before.outlier.removal[
  significant.before.outlier.removal$dispersionOutlier == TRUE,
]


cat(
  "\nGenes passing significance thresholds before",
  "dispersion-outlier removal:",
  nrow(significant.before.outlier.removal),
  "\n"
)

# Genes passing significance thresholds before dispersion-outlier removal: 2646 

cat(
  "Dispersion outliers among these genes:",
  nrow(outliers.among.significant),
  "\n"
)

# Dispersion outliers among these genes: 35 


# 23. Save analysis sample table


write.csv(
  sample.table,
  file = file.path(
    output.dir,
    "RNAseq_DESeq2_sample_table.csv"
  ),
  row.names = FALSE
)



# 24. Save all DESeq2 results


write.csv(
  res.heat.table,
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_all_results.csv"
  ),
  row.names = FALSE
)



# 25. Save final heat-responsive gene set


write.csv(
  heat.sig,
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_DEGs.csv"
  ),
  row.names = FALSE
)



# 26. Save DESeq2 object

# This object contains the fitted model and can be used by
# later scripts without re-running DESeq2.


saveRDS(
  dds.heat,
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_model.rds"
  )
)



# 27. Save model-diagnostic plots


pdf(
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_model_diagnostics.pdf"
  ),
  width = 7,
  height = 6
)


# Dispersion estimates

plotDispEsts(
  dds.heat
)


# MA plot

plotMA(
  res.heat,
  ylim = c(-5, 5)
)


dev.off()



# 28. Save plain-text model summary


capture.output(
  {
    
    cat(
      "RNA-seq acute heat-stress DESeq2 analysis\n"
    )
    
    
    cat(
      "Model formula:\n"
    )
    
    print(
      deseq.formula
    )
    
    cat(
      "\nSamples:",
      ncol(dds.heat),
      "\n"
    )
    
    cat(
      "Genes tested:",
      nrow(dds.heat),
      "\n\n"
    )
    
    cat(
      "Coefficient names:\n"
    )
    
    print(
      resultsNames(dds.heat)
    )
    
    cat(
      "\nFinal coefficient convergence:\n"
    )
    
    print(
      table(
        mcols(dds.heat)$betaConv,
        useNA = "ifany"
      )
    )
    
    cat(
      "\nDispersion outliers:",
      length(high.disp.genes),
      "\n"
    )
    
    cat(
      "\nHeat-stress result summary:\n"
    )
    
    print(
      summary(res.heat)
    )
    
    cat(
      "\nFinal DEG criteria:\n"
    )
    
    cat(
      "Adjusted p-value < 0.01\n"
    )
    
    cat(
      "Absolute log2 fold change > 2\n"
    )
    
    cat(
      "Dispersion outliers excluded\n"
    )
    
    cat(
      "\nFinal heat-responsive genes:",
      number.heat.genes,
      "\n"
    )
    
    cat(
      "Upregulated:",
      number.upregulated,
      "\n"
    )
    
    cat(
      "Downregulated:",
      number.downregulated,
      "\n"
    )
    
  },
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_model_summary.txt"
  )
)



# 29. Final expected-result check


if (nrow(dds.heat) != 19649) {
  
  warning(
    "Expected 19,649 genes in the DESeq2 analysis."
  )
  
}


if (number.heat.genes != 2611) {
  
  warning(
    "Expected 2,611 final heat-responsive genes."
  )
  
}


cat(
  "\nDESeq2 heat-stress analysis complete.\n"
)



# 30. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "RNAseq_heat_stress_sessionInfo.txt"
  )
)
