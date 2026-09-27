# 04_RNAseq_Dixon_validation.R
#
# Purpose:
#   Compare the heat-responsive gene set identified in this
#   study with the Dixon et al. (2020) red co-expression module.
#
#   This script:
#     1. Loads the final heat-responsive gene set from Script 03
#     2. Loads Dixon et al. module assignments
#     3. Tests overlap using a one-sided Fisher's exact test
#     4. Performs PCA on heat-responsive genes
#     5. Performs PCA on Dixon red-module genes
#     6. Calculates mean VST expression of Dixon red-module genes
#        for each sample
#     7. Generates the corresponding expression-validation figure
#
# Inputs:
#   From Script 03:
#     - RNAseq_heat_stress_model.rds
#     - RNAseq_heat_stress_DEGs.csv
#     - RNAseq_DESeq2_sample_table.csv
#
#   External:
#     - wgcna_moduleMembership.tsv
#
# Outputs:
#   - Dixon_red_module_genes.csv
#   - Dixon_heat_overlap_genes.csv
#   - Dixon_heat_overlap_contingency_table.csv
#   - Dixon_heat_overlap_Fisher_test.csv
#   - heat_DEG_PCA_scores.csv
#   - Dixon_red_module_PCA_scores.csv
#   - Dixon_red_module_mean_VST_expression.csv
#   - heat_Dixon_expression_validation.pdf
#   - heat_Dixon_expression_validation.tiff
#   - Dixon_validation_summary.txt
#
# Expected checkpoints:
#   - 634 unique Dixon red-module genes supplied
#   - 611 Dixon red-module genes represented in RNA-seq dataset
#   - 229 genes shared between the Dixon red module and
#     heat-responsive gene set
#   - Fisher's exact test odds ratio approximately 3.15
#   - Dixon red-module PC1 explains approximately 39.1%
#     of expression variation


# 1. Load packages
library(DESeq2)
library(ggplot2)
library(patchwork)


# 2. Define input and output directories


project.dir <- "./Code_and_data/"


# Outputs from Script 03

deseq.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_heat_stress_model.rds"
)

heat.genes.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_heat_stress_DEGs.csv"
)

sample.table.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_DESeq2_sample_table.csv"
)


# Dixon et al. module assignments

dixon.file <- file.path(
  project.dir,
  "./external_resources",
  "wgcna_moduleMembership.tsv"
)


# Output directory

output.dir <- file.path(
  project.dir,
  "./Data/04_Dixon_validation"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# 3. Load DESeq2 model from Script 03

dds.heat <- readRDS(
  deseq.file
)


cat(
  "Genes in final DESeq2 dataset:",
  nrow(dds.heat),
  "\n"
)

# Genes in final DESeq2 dataset: 19649 

cat(
  "Samples in final DESeq2 dataset:",
  ncol(dds.heat),
  "\n"
)

# Samples in final DESeq2 dataset: 58 

# 4. Load final heat-responsive gene set

heat.results <- read.csv(
  heat.genes.file,
  stringsAsFactors = FALSE
)


heat.genes <- unique(
  heat.results$gene
)


cat(
  "Final heat-responsive genes:",
  length(heat.genes),
  "\n"
)

# Final heat-responsive genes: 2611 

# 5. Load sample metadata


sample.table <- read.csv(
  sample.table.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Match metadata order to DESeq2 sample order

sample.table <- sample.table[
  match(
    colnames(dds.heat),
    sample.table$sampleName
  ),
]


if (any(is.na(sample.table$sampleName))) {
  stop(
    "At least one DESeq2 sample was not found in the sample table."
  )
}


if (!all(sample.table$sampleName == colnames(dds.heat))) {
  stop(
    "Sample order does not match between DESeq2 object and sample table."
  )
}


cat(
  "DESeq2 object and sample table are in the same order.\n"
)



# 6. Define heat-stress treatment labels

sample.table$heat_stress_plot <- factor(
  sample.table$heat_stress,
  levels = c(
    "24.5",
    "36"
  ),
  labels = c(
    "24.5°C control",
    "36°C heat stress"
  )
)


cat(
  "\nHeat-stress treatment counts:\n"
)

print(
  table(sample.table$heat_stress_plot)
)
# 24.5°C control 36°C heat stress 
# 29               29 


# PART A: Dixon red-module overlap with heat-responsive genes
# 7. Load Dixon et al. module assignments


dixon.modules <- read.table(
  dixon.file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


cat(
  "\nDixon module-assignment columns:\n"
)

print(
  colnames(dixon.modules)
)
# [1] "gene"          "MMblack"       "MMred"         "MMgrey60"     
# [5] "MMwhite"       "MMdarkred"     "MMgreenyellow" "MMbrown"      
# [9] "MMtan"         "MMlightgreen"  "MMgreen4"      "MMorange"     
# [13] "assignment"   


# 8. Extract Dixon red-module genes

dixon.red.genes <- unique(
  dixon.modules$gene[
    dixon.modules$assignment == "red"
  ]
)


cat(
  "\nUnique Dixon red-module genes supplied:",
  length(dixon.red.genes),
  "\n"
)

# Unique Dixon red-module genes supplied: 634 


# Save complete supplied red-module gene list

write.csv(
  data.frame(
    gene = dixon.red.genes
  ),
  file = file.path(
    output.dir,
    "Dixon_red_module_genes.csv"
  ),
  row.names = FALSE
)



# 9. Define the common gene universe
# The Fisher test should compare only genes that:
#
#   1. were tested in our DESeq2 analysis, AND
#   2. were represented in the Dixon module-assignment dataset.
#
# This ensures that every gene in the contingency table had the
# opportunity to be classified both as heat responsive/not heat
# responsive and red module/not red module.


expression.genes <- rownames(
  dds.heat
)


dixon.assigned.genes <- unique(
  dixon.modules$gene
)


comparison.universe <- intersect(
  expression.genes,
  dixon.assigned.genes
)


cat(
  "\nGenes in common comparison universe:",
  length(comparison.universe),
  "\n"
)

# Genes in common comparison universe: 10408 

# 10. Restrict both gene sets to the common universe


heat.genes.universe <- intersect(
  heat.genes,
  comparison.universe
)


dixon.red.universe <- intersect(
  dixon.red.genes,
  comparison.universe
)


cat(
  "Heat-responsive genes in comparison universe:",
  length(heat.genes.universe),
  "\n"
)

# Heat-responsive genes in comparison universe: 1795 

cat(
  "Dixon red-module genes in comparison universe:",
  length(dixon.red.universe),
  "\n"
)

# Dixon red-module genes in comparison universe: 611 


# 11. Identify genes shared between the two sets

overlap.genes <- intersect(
  heat.genes.universe,
  dixon.red.universe
)


cat(
  "Heat-responsive genes also in Dixon red module:",
  length(overlap.genes),
  "\n"
)
# Heat-responsive genes also in Dixon red module: 229 

# Save overlapping genes

write.csv(
  data.frame(
    gene = overlap.genes
  ),
  file = file.path(
    output.dir,
    "Dixon_heat_overlap_genes.csv"
  ),
  row.names = FALSE
)



# 12. Calculate contingency-table counts


n.universe <- length(
  comparison.universe
)

n.heat <- length(
  heat.genes.universe
)

n.red <- length(
  dixon.red.universe
)

n.overlap <- length(
  overlap.genes
)


# Heat-responsive, but not Dixon red module

n.heat.not.red <- (
  n.heat - n.overlap
)


# Dixon red module, but not heat responsive

n.red.not.heat <- (
  n.red - n.overlap
)


# Neither heat responsive nor Dixon red module

n.neither <- (
  n.universe -
    n.heat -
    n.red +
    n.overlap
)


# 13. Construct 2 x 2 contingency table

overlap.table <- matrix(
  c(
    n.overlap,
    n.heat.not.red,
    n.red.not.heat,
    n.neither
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Heat_response = c(
      "Heat_DEG",
      "Not_heat_DEG"
    ),
    Dixon_module = c(
      "Red_module",
      "Not_red_module"
    )
  )
)


cat(
  "\nDixon red-module overlap contingency table:\n"
)

print(
  overlap.table
)


# Save contingency table
# Dixon_module
# Heat_response  Red_module Not_red_module
# Heat_DEG            229           1566
# Not_heat_DEG        382           8231

write.csv(
  overlap.table,
  file = file.path(
    output.dir,
    "Dixon_heat_overlap_contingency_table.csv"
  ),
  row.names = TRUE
)



# 14. Perform one-sided Fisher's exact test

dixon.fisher <- fisher.test(
  overlap.table,
  alternative = "greater"
)


cat(
  "\nOne-sided Fisher's exact test:\n"
)

print(
  dixon.fisher
)

# Fisher's Exact Test for Count Data
# 
# data:  overlap.table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  2.714217      Inf
# sample estimates:
# odds ratio 
#   3.150424 


# 15. Save Fisher-test summary


dixon.fisher.table <- data.frame(
  comparison_universe = n.universe,
  heat_responsive_genes = n.heat,
  Dixon_red_module_genes = n.red,
  overlapping_genes = n.overlap,
  odds_ratio = unname(
    dixon.fisher$estimate
  ),
  p_value = dixon.fisher$p.value
)


write.csv(
  dixon.fisher.table,
  file = file.path(
    output.dir,
    "Dixon_heat_overlap_Fisher_test.csv"
  ),
  row.names = FALSE
)

# PART B: VST transformation for expression-pattern analyses
# 16. Variance-stabilizing transformation


vst.data <- vst(
  dds.heat,
  blind = TRUE
)


vst.matrix <- assay(
  vst.data
)


cat(
  "\nVST matrix:",
  nrow(vst.matrix),
  "genes x",
  ncol(vst.matrix),
  "samples\n"
)
# VST matrix: 19649 genes x 58 samples


# PART C: PCA of heat-responsive genes
# 17. Identify heat-responsive genes in VST matrix

heat.genes.present <- intersect(
  heat.genes,
  rownames(vst.matrix)
)


cat(
  "\nHeat-responsive genes represented in VST matrix:",
  length(heat.genes.present),
  "\n"
)
# Heat-responsive genes represented in VST matrix: 2611 


# 18. Create heat-responsive-gene expression matrix
# prcomp expects samples in rows and genes in columns.
# The VST matrix needs to be transposed.

heat.expression.matrix <- t(
  vst.matrix[
    heat.genes.present,
    ,
    drop = FALSE
  ]
)



# 19. Perform PCA on heat-responsive genes

heat.pca <- prcomp(
  heat.expression.matrix,
  center = TRUE,
  scale. = TRUE
)


# Percent of variance explained by PC1

heat.pc1.variance <- (
  summary(heat.pca)$importance[2, 1] * 100
)


cat(
  "Heat-responsive-gene PC1 variance explained:",
  round(heat.pc1.variance, 1),
  "%\n"
)

# Heat-responsive-gene PC1 variance explained: 55.1 %

# 20. Save heat-responsive-gene PC1 scores


heat.pca.scores <- data.frame(
  sample = rownames(heat.pca$x),
  PC1 = heat.pca$x[, 1],
  heat_stress = sample.table$heat_stress_plot,
  stringsAsFactors = FALSE
)


write.csv(
  heat.pca.scores,
  file = file.path(
    output.dir,
    "heat_DEG_PCA_scores.csv"
  ),
  row.names = FALSE
)



# PART D: PCA of Dixon red-module genes
# 21. Identify Dixon red-module genes in VST matrix

dixon.red.genes.present <- intersect(
  dixon.red.genes,
  rownames(vst.matrix)
)


cat(
  "\nDixon red-module genes represented in VST matrix:",
  length(dixon.red.genes.present),
  "\n"
)


# 22. Create Dixon red-module expression matrix

dixon.expression.matrix <- t(
  vst.matrix[
    dixon.red.genes.present,
    ,
    drop = FALSE
  ]
)



# 23. Perform PCA on Dixon red-module genes


dixon.pca <- prcomp(
  dixon.expression.matrix,
  center = TRUE,
  scale. = TRUE
)


# Percent of variance explained by PC1

dixon.pc1.variance <- (
  summary(dixon.pca)$importance[2, 1] * 100
)


cat(
  "Dixon red-module PC1 variance explained:",
  round(dixon.pc1.variance, 1),
  "%\n"
)
# Dixon red-module PC1 variance explained: 39.1 %


# 24. Save Dixon red-module PC1 scores

dixon.pca.scores <- data.frame(
  sample = rownames(dixon.pca$x),
  PC1 = dixon.pca$x[, 1],
  heat_stress = sample.table$heat_stress_plot,
  stringsAsFactors = FALSE
)


write.csv(
  dixon.pca.scores,
  file = file.path(
    output.dir,
    "Dixon_red_module_PCA_scores.csv"
  ),
  row.names = FALSE
)

# PART E: Mean Dixon red-module VST expression

# 25. Calculate mean red-module expression for each sample
# For each sample, calculate the mean VST-normalized expression
# across all Dixon red-module genes represented in the
# expression matrix.


dixon.mean.vst <- colMeans(
  vst.matrix[
    dixon.red.genes.present,
    ,
    drop = FALSE
  ]
)


dixon.mean.vst.data <- data.frame(
  sample = names(dixon.mean.vst),
  mean_vst_expression = as.numeric(dixon.mean.vst),
  heat_stress = sample.table$heat_stress_plot,
  stringsAsFactors = FALSE
)


write.csv(
  dixon.mean.vst.data,
  file = file.path(
    output.dir,
    "Dixon_red_module_mean_VST_expression.csv"
  ),
  row.names = FALSE
)

# PART F: Create expression-validation figure
# 26. Combine PCA scores for plotting

heat.pca.plot.data <- data.frame(
  sample = heat.pca.scores$sample,
  PC1 = heat.pca.scores$PC1,
  heat_stress = heat.pca.scores$heat_stress,
  module = "Heat-responsive genes",
  stringsAsFactors = FALSE
)


dixon.pca.plot.data <- data.frame(
  sample = dixon.pca.scores$sample,
  PC1 = dixon.pca.scores$PC1,
  heat_stress = dixon.pca.scores$heat_stress,
  module = "Dixon red module",
  stringsAsFactors = FALSE
)


pca.plot.data <- rbind(
  heat.pca.plot.data,
  dixon.pca.plot.data
)



# 27. Add variance explained to PCA facet labels

pca.plot.data$facet_label[
  pca.plot.data$module == "Heat-responsive genes"
] <- paste0(
  "Heat-responsive genes",
  "\nPC1 = ",
  round(heat.pc1.variance, 1),
  "% variance explained"
)


pca.plot.data$facet_label[
  pca.plot.data$module == "Dixon red module"
] <- paste0(
  "Dixon red module",
  "\nPC1 = ",
  round(dixon.pc1.variance, 1),
  "% variance explained"
)


# Control facet order

pca.plot.data$facet_label <- factor(
  pca.plot.data$facet_label,
  levels = c(
    paste0(
      "Heat-responsive genes",
      "\nPC1 = ",
      round(heat.pc1.variance, 1),
      "% variance explained"
    ),
    paste0(
      "Dixon red module",
      "\nPC1 = ",
      round(dixon.pc1.variance, 1),
      "% variance explained"
    )
  )
)


# 28. Define treatment colors


treatment.colors <- c(
  "24.5°C control" = "#bfd7ea",
  "36°C heat stress" = "#ff5a5f"
)



# 29. Plot PCA results

pca.module.plot <- ggplot(
  pca.plot.data,
  aes(
    x = heat_stress,
    y = PC1,
    fill = heat_stress,
    color = heat_stress
  )
) +
  geom_boxplot(
    alpha = 0.6,
    width = 0.6,
    outlier.shape = NA,
    linewidth = 1
  ) +
  geom_jitter(
    shape = 21,
    size = 3,
    alpha = 0.75,
    width = 0.12,
    height = 0,
    stroke = 0.8
  ) +
  facet_wrap(
    ~ facet_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = treatment.colors,
    name = "Treatment"
  ) +
  scale_color_manual(
    values = treatment.colors,
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = "PC1 score"
  ) +
  theme_bw(
    base_size = 20
  ) +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      color = "grey40",
      linewidth = 0.8
    ),
    strip.text = element_text(
      size = 18,
      face = "bold",
      margin = margin(
        8,
        8,
        8,
        8
      )
    ),
    axis.text.x = element_text(
      size = 17,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 16
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )


pca.module.plot



# 30. Add facet label to mean-expression panel

dixon.mean.vst.data$facet_label <- paste0(
  "Dixon red module",
  "\nMean VST expression across ",
  length(dixon.red.genes.present),
  " genes"
)



# 31. Plot mean Dixon red-module expression


dixon.mean.vst.plot <- ggplot(
  dixon.mean.vst.data,
  aes(
    x = heat_stress,
    y = mean_vst_expression,
    fill = heat_stress,
    color = heat_stress
  )
) +
  geom_boxplot(
    alpha = 0.6,
    width = 0.6,
    outlier.shape = NA,
    linewidth = 1
  ) +
  geom_jitter(
    shape = 21,
    size = 3,
    alpha = 0.75,
    width = 0.12,
    height = 0,
    stroke = 0.8
  ) +
  facet_wrap(
    ~ facet_label,
    nrow = 1
  ) +
  scale_fill_manual(
    values = treatment.colors,
    name = "Treatment"
  ) +
  scale_color_manual(
    values = treatment.colors,
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = "Mean VST expression"
  ) +
  theme_bw(
    base_size = 20
  ) +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      color = "grey40",
      linewidth = 0.8
    ),
    strip.text = element_text(
      size = 18,
      face = "bold",
      margin = margin(
        8,
        8,
        8,
        8
      )
    ),
    axis.text.x = element_text(
      size = 17,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 16
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )


dixon.mean.vst.plot



# 32. Combine PCA and mean-expression panels

expression.validation.plot <- (
  pca.module.plot |
    dixon.mean.vst.plot
) +
  plot_layout(
    widths = c(
      2,
      1
    )
  )


expression.validation.plot



# 33. Save combined figure

ggsave(
  filename = file.path(
    output.dir,
    "heat_Dixon_expression_validation.pdf"
  ),
  plot = expression.validation.plot,
  width = 18,
  height = 7,
  units = "in"
)


ggsave(
  filename = file.path(
    output.dir,
    "heat_Dixon_expression_validation.tiff"
  ),
  plot = expression.validation.plot,
  width = 18,
  height = 7,
  units = "in",
  dpi = 300,
  compression = "lzw"
)



# PART G: Final checks and summary
# 34. Print expected checkpoints


cat(
  "Unique Dixon red-module genes supplied:",
  length(dixon.red.genes),
  "\n"
)
# Unique Dixon red-module genes supplied: 634 

cat(
  "Dixon red-module genes represented in expression dataset:",
  length(dixon.red.genes.present),
  "\n"
)
# Dixon red-module genes represented in expression dataset: 611 

cat(
  "Heat-responsive genes overlapping Dixon red module:",
  length(overlap.genes),
  "\n"
)
# Heat-responsive genes overlapping Dixon red module: 229 

cat(
  "Fisher's exact-test odds ratio:",
  unname(dixon.fisher$estimate),
  "\n"
)
# Fisher's exact-test odds ratio: 3.150424 

cat(
  "Fisher's exact-test p-value:",
  dixon.fisher$p.value,
  "\n"
)
# Fisher's exact-test p-value: 4.084871e-35 

cat(
  "Dixon red-module PC1 variance explained:",
  round(dixon.pc1.variance, 1),
  "%\n"
)

# Dixon red-module PC1 variance explained: 39.1 %



# 35. Save plain-text analysis summary


capture.output(
  {
    
    cat(
      "Dixon red-module validation\n"
    )
    

    
    cat(
      "Unique Dixon red-module genes supplied:",
      length(dixon.red.genes),
      "\n"
    )
    
    cat(
      "Dixon red-module genes represented in RNA-seq dataset:",
      length(dixon.red.genes.present),
      "\n"
    )
    
    cat(
      "Genes in Fisher-test comparison universe:",
      length(comparison.universe),
      "\n"
    )
    
    cat(
      "Heat-responsive genes in comparison universe:",
      length(heat.genes.universe),
      "\n"
    )
    
    cat(
      "Dixon red-module genes in comparison universe:",
      length(dixon.red.universe),
      "\n"
    )
    
    cat(
      "Heat-responsive / Dixon red-module overlap:",
      length(overlap.genes),
      "\n\n"
    )
    
    cat(
      "One-sided Fisher's exact test:\n"
    )
    
    print(
      dixon.fisher
    )
    
    cat(
      "\nHeat-responsive-gene PC1 variance explained:",
      round(heat.pc1.variance, 1),
      "%\n"
    )
    
    cat(
      "Dixon red-module PC1 variance explained:",
      round(dixon.pc1.variance, 1),
      "%\n"
    )
    
  },
  file = file.path(
    output.dir,
    "Dixon_validation_summary.txt"
  )
)



# 37. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "Dixon_validation_sessionInfo.txt"
  )
)
