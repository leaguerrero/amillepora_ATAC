# 09_ATAC_RNA_integration.R
#
# Purpose:
#   Examine broad associations between baseline promoter
#   accessibility and gene-expression characteristics using
#   independently generated ATAC-seq and RNA-seq datasets from
#   Acropora millepora.
#
# Important interpretation:
#   The ATAC-seq and RNA-seq datasets were generated from
#   independently sourced fragments, and genotype correspondence
#   between the experiments is unknown.
#
#   Therefore, these analyses describe broad cross-dataset
#   associations. 
#
#
# This script:
#
#   1. Loads ATAC-seq consensus-peak raw counts and annotations
#   2. Calculates CPM from the four ATAC-seq sample columns only
#   3. Calculates mean promoter accessibility across those four samples
#   4. Loads the final RNA-seq DESeq2 model
#   5. Calculates mean normalized expression across 58 samples
#   6. Integrates genes represented in both datasets
#   7. Retains the promoter peak closest to the TSS for each gene
#
#   Figure 3A:
#     promoter accessibility by mean-expression quintile
#
#   Figure 3B:
#     promoter accessibility by expression-variability quintile
#
#   Promoter enrichment:
#     - heat-responsive genes
#     - Dixon et al. red-module genes
#
#   Figure 4:
#     promoter accessibility vs:
#       - mean normalized expression
#       - magnitude of heat-induced expression change
#
#
# Inputs:
#
#   Script 03:
#     RNAseq_heat_stress_model.rds
#     RNAseq_heat_stress_DEGs.csv
#
#   Script 06:
#     ATAC_consensus_peak_counts.csv
#     ATAC_annotated_consensus_peaks.csv
#
#   Dixon et al. module membership:
#     wgcna_moduleMembership.tsv
#
#
# Expected manuscript checkpoints:
#
#   19,649 genes in RNA-seq expression universe
#   2,611 heat-responsive genes
#   1,735 expressed genes with accessible promoters
#
#   Heat-responsive genes:
#     243 with accessible promoters
#     Fisher OR ~ 1.07
#     p ~ 0.355
#
#   Dixon red module:
#     611 genes represented in expression universe
#     52 with accessible promoters
#     Fisher OR ~ 0.96
#     p ~ 0.828
#
#   Figure 4 correlation values are recalculated below from the
#   corrected ATAC-seq accessibility metric.

# 1. Load packages

library(DESeq2)
library(edgeR)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(GenomicRanges)

# 2. Define project directory
# Run this script from the Code_and_data directory.
project.dir <- "."

# 3. Define input files

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


atac.counts.file <- file.path(
  project.dir,
  "Data/06_consensus_peaks",
  "ATAC_consensus_peak_counts.csv"
)

atac.annotations.file <- file.path(
  project.dir,
  "Data/06_consensus_peaks",
  "ATAC_annotated_consensus_peaks.csv"
)

atac.consensus.anno.file <- file.path(
  project.dir,
  "Data/06_consensus_peaks",
  "ATAC_consensus_peaks_with_regions.csv"
)


dixon.file <- file.path(
  project.dir,
  "./external_resources",
  "wgcna_moduleMembership.tsv"
)

snp.file <- file.path(
  project.dir,
  "Data/09_RNA_integration",
  "out.snpden"
)

# 4. Define output directory

output.dir <- file.path(
  project.dir,
  "./Data/09_RNA_integration"
)


dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# PART A: Load RNA-seq data

# 5. Load final heat-stress DESeq2 object
dds.heat <- readRDS(
  deseq.file
)


cat(
  "RNA-seq genes in DESeq2 object:",
  nrow(dds.heat),
  "\n"
)



# 6. Load final heat-responsive gene set


heat.results <- read.csv(
  heat.genes.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)




heat.results$gene <- as.character(
  heat.results$gene
)


heat.results <- heat.results[
  !duplicated(heat.results$gene),
]

# 7. Extract DESeq2-normalized expression counts

normalized.counts <- counts(
  dds.heat,
  normalized = TRUE
)


dim(normalized.counts)

# [1] 19649    58


# 8. Calculate mean normalized expression
# DESeq2 baseMean is the mean of normalized counts across
# samples.

# We calculate it directly here so the value used in the
# integration analysis is transparent.


mean.expression <- rowMeans(
  normalized.counts
)


expression.data <- data.frame(
  gene = rownames(normalized.counts),
  baseMean = mean.expression,
  stringsAsFactors = FALSE
)

summary(expression.data$baseMean)

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.06     0.71     3.61    38.45    17.33 51167.57 


# PART B: Calculate corrected ATAC-seq accessibility

# 9. Load raw counts for the 11,231 consensus peaks
#
# IMPORTANT:
# ATAC_consensus_peak_counts.csv contains genomic-coordinate
# metadata plus four ATAC-seq sample-count columns. Only the
# four biological-sample count columns are normalized.
#
# Accessibility is defined as the mean CPM across ATAC-seq
# samples 1-4 for each consensus peak.

atac.counts.table <- read.csv(
  atac.counts.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

nrow(atac.counts.table)
# [1] 11231

# 10. Confirm required count-table columns

atac.sample.columns <- c(
  "1",
  "2",
  "3",
  "4"
)

required.count.columns <- c(
  "peak_id",
  atac.sample.columns
)


# 11. Extract ONLY the four biological-sample count columns

atac.raw.counts <- as.matrix(
  atac.counts.table[
    ,
    atac.sample.columns,
    drop = FALSE
  ]
)

storage.mode(atac.raw.counts) <- "numeric"

dim(atac.raw.counts)
# Expected: 11231 x 4
# [1] 11231     4

# 12. CPM-normalize the four ATAC-seq samples

atac.cpm <- edgeR::cpm(
  atac.raw.counts
)

dim(atac.cpm)

# Expected: 11231 x 4


# 13. Calculate mean normalized CPM for each consensus peak

mean_normalized_CPM <- rowMeans(
  atac.cpm
)

summary(mean_normalized_CPM)

# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.52     14.72     32.47     89.04     75.75 161996.90 

# 14. Build peak-level accessibility table

peak.accessibility <- data.frame(
  peak_id = as.character(
    atac.counts.table$peak_id
  ),
  mean_normalized_CPM = mean_normalized_CPM,
  stringsAsFactors = FALSE
)


# 15. Load ChIPseeker annotations
#
# The legacy avgAccess column is explicitly discarded. Corrected
# accessibility is joined by peak_id instead of row position.

annotated.peaks <- read.csv(
  atac.annotations.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


annotated.peaks <- annotated.peaks %>%
  dplyr::select(
    -any_of(
      c(
        "avgAccess",
        "avgAccessibility",
        "mean_CPM",
        "mean_normalized_CPM"
      )
    )
  ) %>%
  left_join(
    peak.accessibility,
    by = "peak_id"
  )

nrow(annotated.peaks)
# [1] 11096


# 16. Confirm required ATAC columns

required.atac.columns <- c(
  "peak_id",
  "gene",
  "annotation",
  "distanceToTSS",
  "mean_normalized_CPM"
)

missing.atac.columns <- setdiff(
  required.atac.columns,
  colnames(annotated.peaks)
)

# 17. Ensure gene identifiers are character values

annotated.peaks$gene <- as.character(
  annotated.peaks$gene
)


# PART C: Integrate ATAC-seq and RNA-seq data
# 18. Merge expression data with ATAC-seq peak annotations

all.genes.and.peaks <- expression.data %>%
  inner_join(
    annotated.peaks,
    by = "gene"
  )

nrow(all.genes.and.peaks)
# [1] 4985

# Rows with both RNA-seq and ATAC-seq information: 4985 
n_distinct(all.genes.and.peaks$gene)
# [1] 3638

# Unique genes with both RNA-seq and ATAC-seq information: 3638 

# 19. Flag promoter-associated peaks

all.genes.and.peaks <- all.genes.and.peaks %>%
  mutate(
    Promoter = grepl(
      "promoter",
      annotation,
      ignore.case = TRUE
    )
  )

sum(all.genes.and.peaks$Promoter,
    na.rm = TRUE)
# [1] 1919
# Promoter-associated peak rows in integrated dataset: 1919 

# PART D: Select one promoter peak per gene
# 20. Retain promoter peak closest to TSS
# For genes associated with multiple promoter-region consensus
# peaks, retain the peak with the smallest absolute distance
# to the transcription start site.
#
# This is the gene-level promoter-accessibility value used
# throughout the remaining integration analyses.


closest.promoter.per.gene <- all.genes.and.peaks %>%
  filter(
    Promoter,
    !is.na(gene),
    is.finite(distanceToTSS),
    is.finite(mean_normalized_CPM),
    is.finite(baseMean)
  ) %>%
  mutate(
    abs.distanceToTSS = abs(
      distanceToTSS
    )
  ) %>%
  group_by(
    gene
  ) %>%
  slice_min(
    order_by = abs.distanceToTSS,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()



# 21. Confirm one observation per gene
sum(all.genes.and.peaks$Promoter, na.rm = TRUE)
# [1] 1919

# Promoter-associated peak rows before selecting closest peak: 1919 

nrow(closest.promoter.per.gene)
# [1] 1735


sum(duplicated(closest.promoter.per.gene$gene))
# [1] 0


# Duplicated genes remaining: 0 

# Current manuscript reports 1,735 genes


# 22. Save integrated closest-promoter table


write.csv(
  closest.promoter.per.gene,
  file = file.path(
    output.dir,
    "ATAC_RNA_closest_promoter_per_gene.csv"
  ),
  row.names = FALSE
)


# PART E: Figure 3B
# Mean expression and promoter accessibility

# 23. Define gene-level expression quintiles
# expression quintiles are calculated among genes represented
# in the integrated ATAC/RNA dataset.


expression.quintile.data <- all.genes.and.peaks %>%
  filter(
    !is.na(gene),
    is.finite(baseMean)
  ) %>%
  distinct(
    gene,
    .keep_all = TRUE
  )



# 24. Calculate expression-quintile boundaries
expression.breaks <- quantile(
  expression.quintile.data$baseMean,
  probs = seq(
    0,
    1,
    by = 0.2
  ),
  na.rm = TRUE
)


expression.breaks
# 0%          20%          40%          60%          80%         100% 
#   7.918801e-02 3.704202e-01 1.058396e+00 3.788327e+00 1.869896e+01 5.116757e+04 


# 25. Assign expression quintiles
expression.quintile.data <- expression.quintile.data %>%
  mutate(
    expression.quintile = cut(
      baseMean,
      breaks = expression.breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
  ) %>%
  dplyr::select(
    gene,
    expression.quintile
  )

table(expression.quintile.data$expression.quintile)

# 1   2   3   4   5 
# 728 727 728 727 728 

# 26. Join expression quintiles to closest-promoter data

promoter.expression.data <- closest.promoter.per.gene %>%
  left_join(
    expression.quintile.data,
    by = "gene"
  ) %>%
  filter(
    expression.quintile %in% c(
      1,
      5
    )
  ) %>%
  mutate(
    expression.category = factor(
      if_else(
        expression.quintile == 1,
        "Low expression",
        "High expression"
      ),
      levels = c(
        "Low expression",
        "High expression"
      )
    )
  )

table(promoter.expression.data$expression.category)

# Low expression High expression 
# 360             362 

# 27. Build Figure 3B


promoter.expression.plot <- ggplot(
  promoter.expression.data,
  aes(
    x = distanceToTSS,
    y = mean_normalized_CPM,
    color = expression.category,
    fill = expression.category
  )
) +
  geom_smooth(
    se = TRUE,
    linewidth = 1.2
  ) +
  labs(
    x = "Distance to TSS (bp)",
    y = "Accessibility\n(mean normalized CPM)"
  ) +
  scale_color_manual(
    name = "Expression category",
    values = c(
      "Low expression" = "#16db93",
      "High expression" = "#f29e4c"
    )
  ) +
  scale_fill_manual(
    name = "Expression category",
    values = c(
      "Low expression" = "#C7F9E7",
      "High expression" = "#FBE0C6"
    )
  ) +
  theme_pubr(
    base_size = 20
  ) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "plain",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 20,
      face = "plain"
    ),
    axis.text = element_text(
      size = 17
    ),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.16),
    legend.title = element_text(
      size = 16,
      face = "plain"
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.key.height = unit(0.55, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.background = element_rect(
      fill = "white",
      color = NA
    )
  ) +
  coord_cartesian(
    xlim = c(
      -1000,
      1000
    )
  )


promoter.expression.plot

# PART F: Figure 3C
# Expression variability and promoter accessibility

# 28. Calculate gene-level expression standard deviation

expression.standard.deviation <- apply(
  normalized.counts,
  1,
  sd,
  na.rm = TRUE
)



# 29. Build expression-variability table


expression.variation.data <- data.frame(
  gene = rownames(normalized.counts),
  baseMean = mean.expression,
  standard.deviation = expression.standard.deviation,
  stringsAsFactors = FALSE
)



# 30. Exclude extreme-SD genes


# Existing analysis excludes genes with:
#
#     standard deviation > 20,000

expression.variation.data <- expression.variation.data %>%
  filter(
    is.finite(standard.deviation),
    standard.deviation <= 20000,
    is.finite(baseMean),
    baseMean > 0
  )

nrow(expression.variation.data)
# [1] 19648



# 31. Calculate coefficient of variation

expression.variation.data <- expression.variation.data %>%
  mutate(
    coefficient.of.variation =
      standard.deviation / baseMean
  ) %>%
  filter(
    is.finite(
      coefficient.of.variation
    )
  )



# 32. Join expression variation to closest-promoter genes

promoter.variation.data <- closest.promoter.per.gene %>%
  dplyr::select(
    -baseMean
  ) %>%
  inner_join(
    expression.variation.data,
    by = "gene"
  )


nrow(promoter.variation.data)
# [1] 1734

# Promoter-accessible genes retained for CV analysis: 1734 

# 33. Calculate coefficient-of-variation quintiles

cv.breaks <- quantile(
  promoter.variation.data$coefficient.of.variation,
  probs = seq(
    0,
    1,
    by = 0.2
  ),
  na.rm = TRUE
)

cv.breaks

# 0%       20%       40%       60%       80%      100% 
#   0.2181854 0.7929374 1.1959414 1.7576145 2.4330158 6.0590039 


# 34. Assign CV quintiles
promoter.variation.data <- promoter.variation.data %>%
  mutate(
    variability.quintile = cut(
      coefficient.of.variation,
      breaks = cv.breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
  )


table(promoter.variation.data$variability.quintile)

# 1   2   3   4   5 
# 347 347 346 347 347 

# 35. Retain lowest and highest variability quintiles

promoter.variation.extremes <- promoter.variation.data %>%
  filter(
    variability.quintile %in% c(
      1,
      5
    )
  ) %>%
  mutate(
    expression.variation = factor(
      if_else(
        variability.quintile == 1,
        "Low variability",
        "High variability"
      ),
      levels = c(
        "Low variability",
        "High variability"
      )
    )
  )


table(promoter.variation.extremes$expression.variation)

# Low variability High variability 
# 347              347 

# 31. Build Figure 3C


promoter.variation.plot <- ggplot(
  promoter.variation.extremes,
  aes(
    x = distanceToTSS,
    y = mean_normalized_CPM,
    color = expression.variation,
    fill = expression.variation
  )
) +
  geom_smooth(
    se = TRUE,
    linewidth = 1.2
  ) +
  labs(
    x = "Distance to TSS (bp)",
    y = "Accessibility\n(mean normalized CPM)"
  ) +
  scale_color_manual(
    name = "Expression variation",
    values = c(
      "Low variability" = "#8591FF",
      "High variability" = "#97d9e4"
    )
  ) +
  scale_fill_manual(
    name = "Expression variation",
    values = c(
      "Low variability" = "#aab2ff",
      "High variability" = "#CEEDF2"
    )
  ) +
  theme_pubr(
    base_size = 20
  ) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "plain",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 20,
      face = "plain"
    ),
    axis.text = element_text(
      size = 17
    ),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.16),
    legend.title = element_text(
      size = 16,
      face = "plain"
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.key.height = unit(0.55, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.background = element_rect(
      fill = "white",
      color = NA
    )
  ) +
  coord_cartesian(
    xlim = c(
      -1000,
      1000
    )
  )


promoter.variation.plot


# PART G: SNP Density plot

##read in data
# atac <- read.csv("~/Desktop/ATAC_consensus_peaks_with_regions.csv")
atac.consensus.anno <- read.csv(atac.consensus.anno.file)
atac.consensus.anno <- separate(atac.consensus.anno,col="annotation",into=c("type","more"),extra="merge")
atac.consensus.anno.filt <- atac.consensus.anno %>% 
  filter(type%in%c("Distal","Intron","Promoter"))

snp <- read.delim(snp.file)
snp$BIN_END <- snp$BIN_START+999
snp$name <- paste(snp$CHROM,snp$BIN_START)

##make granges objects
snp.gr <- GRanges(seqnames=snp$CHROM,
                  ranges=IRanges(start=snp$BIN_START,end=snp$BIN_END))
atac.gr <- GRanges(seqnames=atac.consensus.anno.filt$chromosome,
                   ranges=IRanges(start=atac.consensus.anno.filt$peak_center,
                                  end=atac.consensus.anno.filt$peak_center))

##find nearest peak
nearest_peak <- distanceToNearest(snp.gr,atac.gr)
bin_mid <- start(snp.gr[queryHits(nearest_peak)]) + 500
peak_pos <- start(atac.gr)[subjectHits(nearest_peak)]

##Signed dist
signed_dist <- peak_pos-bin_mid

snp$dist <- NA
snp$dist[queryHits(nearest_peak)] <- signed_dist
snp$type <- NA
snp$type[queryHits(nearest_peak)] <- atac.consensus.anno.filt$type[subjectHits(nearest_peak)]
head(snp)


snp.plot <- ggplot(
  snp %>%
    filter(abs(dist) < 100000) %>%
    mutate(
      type = factor(
        type,
        levels = c(
          "Distal",
          "Intron",
          "Promoter"
        )
      )
    ),
  aes(
    x = dist,
    y = SNP_COUNT,
    color = type,
    fill = type
  )
) +
  geom_smooth(
    se = TRUE,
    linewidth = 1.2
  ) +
  labs(
    x = "Distance to peak (bp)",
    y = "SNP density\n(SNPs / kb)"
  ) +
  scale_color_manual(
    name = "Region type",
    values = c(
      "Distal" = "#0e0e52",
      "Intron" = "#679436",
      "Promoter" = "#ff5400"
    )
  ) +
  scale_fill_manual(
    name = "Region type",
    values = c(
      "Distal" = "#449dd1",
      "Intron" = "#a5be00",
      "Promoter" = "#ff9e00"
    )
  ) +
  theme_pubr(
    legend = "bottom",
    base_size = 20
  ) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "plain",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 20,
      face = "plain"
    ),
    axis.text = element_text(
      size = 17
    ),
    legend.title = element_text(
      size = 16,
      face = "plain"
    ),
    legend.text = element_text(
      size = 15
    )
  ) +
  coord_cartesian(
    xlim = c(
      -100000,
      100000
    )
  )

snp.plot

# PART H: Save Figure 3
# 32. Combine Figure 3 panels


figure3 <- snp.plot + (promoter.expression.plot /
  promoter.variation.plot) +
  plot_annotation(
    tag_levels = "A"
  )


figure3



# 33. Save Figure 3

ggsave(
  filename = file.path(
    output.dir,
    "Figure3_ATAC_expression_associations.pdf"
  ),
  plot = figure3,
  width = 17,
  height = 13
)


ggsave(
  filename = file.path(
    output.dir,
    "Figure3_ATAC_expression_associations.tiff"
  ),
  plot = figure3,
  width = 15,
  height = 7,
  dpi = 300,
  compression = "lzw"
)



# PART I: Heat-responsive gene promoter enrichment
# 34. Define RNA-seq expression universe

expression.universe <- unique(
  rownames(dds.heat)
)

length(expression.universe)
# [1] 19649

# Expression-analysis universe: 19649 


# 35. Define heat-responsive genes in expression universe
heat.genes <- intersect(
  heat.results$gene,
  expression.universe
)


length(heat.genes)
# [1] 2611

# Heat-responsive genes in expression universe: 2611 

# 36. Define accessible-promoter genes in expression universe

accessible.promoter.genes <- annotated.peaks %>%
  filter(
    grepl(
      "promoter",
      annotation,
      ignore.case = TRUE
    ),
    !is.na(gene),
    gene != ""
  ) %>%
  pull(
    gene
  ) %>%
  unique() %>%
  intersect(
    expression.universe
  )

length(accessible.promoter.genes)
# [1] 1735

# Genes with accessible promoters in expression universe: 1735 


# 37. Define non-heat-responsive genes

nonheat.genes <- setdiff(
  expression.universe,
  heat.genes
)


# 38. Count heat-responsive genes with accessible promoters

heat.accessible <- length(
  intersect(
    heat.genes,
    accessible.promoter.genes
  )
)


# 39. Count heat-responsive genes without accessible promoters


heat.not.accessible <- length(
  setdiff(
    heat.genes,
    accessible.promoter.genes
  )
)



# 40. Count non-heat-responsive genes with accessible promoters

nonheat.accessible <- length(
  intersect(
    nonheat.genes,
    accessible.promoter.genes
  )
)



# 41. Count non-heat-responsive genes without accessible promoters

nonheat.not.accessible <- length(
  setdiff(
    nonheat.genes,
    accessible.promoter.genes
  )
)



# 42. Build heat-response contingency table

heat.promoter.table <- matrix(
  c(
    heat.accessible,
    heat.not.accessible,
    nonheat.accessible,
    nonheat.not.accessible
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Heat_response = c(
      "Heat_DEG",
      "Not_heat_DEG"
    ),
    Promoter = c(
      "Accessible",
      "Not_accessible"
    )
  )
)


heat.promoter.table

# Promoter
# Heat_response  Accessible Not_accessible
# Heat_DEG            243           2368
# Not_heat_DEG       1492          15546


# 43. Run two-sided Fisher's exact test

heat.promoter.fisher <- fisher.test(
  heat.promoter.table,
  alternative = "two.sided"
)


heat.promoter.fisher


# Fisher's Exact Test for Count Data
# 
# data:  heat.promoter.table
# p-value = 0.3546
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9235019 1.2340791
# sample estimates:
# odds ratio 
#   1.069248 


# 44. Save heat-response Fisher result


heat.promoter.summary <- data.frame(
  expression_universe = length(expression.universe),
  heat_responsive_genes = length(heat.genes),
  heat_responsive_accessible = heat.accessible,
  nonheat_genes = length(nonheat.genes),
  nonheat_accessible = nonheat.accessible,
  odds_ratio = unname(
    heat.promoter.fisher$estimate
  ),
  ci_lower = heat.promoter.fisher$conf.int[1],
  ci_upper = heat.promoter.fisher$conf.int[2],
  p_value = heat.promoter.fisher$p.value
)


write.csv(
  heat.promoter.summary,
  file = file.path(
    output.dir,
    "ATAC_heat_response_promoter_Fisher.csv"
  ),
  row.names = FALSE
)



# PART J: Dixon red-module promoter enrichment
# 45. Load Dixon et al. module membership

dixon.modules <- read.table(
  dixon.file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# 46. Extract Dixon red-module genes


dixon.red <- dixon.modules %>%
  filter(
    assignment == "red"
  ) %>%
  pull(
    gene
  ) %>%
  unique()


length(dixon.red)
# [1] 634

# Total unique Dixon red-module genes supplied: 634 


# 47. Restrict Dixon module to RNA-seq expression universe

dixon.red.universe <- intersect(
  dixon.red,
  expression.universe
)

length(dixon.red.universe)
# [1] 611

# Dixon red-module genes represented in expression universe: 611 

# 48. Define genes outside Dixon red module

not.dixon.red <- setdiff(
  expression.universe,
  dixon.red.universe
)


# 49. Count red-module genes with accessible promoters

dixon.accessible <- length(
  intersect(
    dixon.red.universe,
    accessible.promoter.genes
  )
)


# 50. Count red-module genes without accessible promoters

dixon.not.accessible <- length(
  setdiff(
    dixon.red.universe,
    accessible.promoter.genes
  )
)



# 51. Count non-red genes with accessible promoters

nondixon.accessible <- length(
  intersect(
    not.dixon.red,
    accessible.promoter.genes
  )
)


# 52. Count non-red genes without accessible promoters

nondixon.not.accessible <- length(
  setdiff(
    not.dixon.red,
    accessible.promoter.genes
  )
)

# 53. Build Dixon contingency table

dixon.promoter.table <- matrix(
  c(
    dixon.accessible,
    dixon.not.accessible,
    nondixon.accessible,
    nondixon.not.accessible
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Dixon_module = c(
      "Red_module",
      "Not_red_module"
    ),
    Promoter = c(
      "Accessible",
      "Not_accessible"
    )
  )
)

dixon.promoter.table


# Promoter
# Dixon_module     Accessible Not_accessible
# Red_module             52            559
# Not_red_module       1683          17355


# 54. Run two-sided Fisher's exact test

dixon.promoter.fisher <- fisher.test(
  dixon.promoter.table,
  alternative = "two.sided"
)


dixon.promoter.fisher


# Fisher's Exact Test for Count Data
# 
# data:  dixon.promoter.table
# p-value = 0.8282
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7044794 1.2823164
# sample estimates:
# odds ratio 
#  0.9592505 

# 55. Save Dixon Fisher result

dixon.promoter.summary <- data.frame(
  expression_universe = length(expression.universe),
  red_module_genes = length(dixon.red.universe),
  red_module_accessible = dixon.accessible,
  nonred_genes = length(not.dixon.red),
  nonred_accessible = nondixon.accessible,
  odds_ratio = unname(
    dixon.promoter.fisher$estimate
  ),
  ci_lower = dixon.promoter.fisher$conf.int[1],
  ci_upper = dixon.promoter.fisher$conf.int[2],
  p_value = dixon.promoter.fisher$p.value
)


write.csv(
  dixon.promoter.summary,
  file = file.path(
    output.dir,
    "ATAC_Dixon_red_module_promoter_Fisher.csv"
  ),
  row.names = FALSE
)



# PART K: Heat-responsive genes with accessible promoters


# 56. Prepare final heat-responsive promoter dataset
#
# The same closest-promoter-per-gene definition used above is
# used here.



heat.lfc.data <- heat.results %>%
  dplyr::select(
    gene,
    log2FoldChange
  )


heat.promoter.data <- closest.promoter.per.gene %>%
  inner_join(
    heat.lfc.data,
    by = "gene"
  ) %>%
  filter(
    is.finite(mean_normalized_CPM),
    mean_normalized_CPM > 0,
    is.finite(baseMean),
    baseMean > 0,
    is.finite(log2FoldChange)
  )


nrow(heat.promoter.data)
# [1] 243

# Heat-responsive genes with usable accessible-promoter data: 243 


# 57. Save Figure 4 input data


write.csv(
  heat.promoter.data,
  file = file.path(
    output.dir,
    "ATAC_heat_response_accessible_promoter_data.csv"
  ),
  row.names = FALSE
)



# PART L: Spearman correlations

# 58. Promoter accessibility vs mean expression


cor.mean.expression <- cor.test(
  log10(
    heat.promoter.data$mean_normalized_CPM
  ),
  log10(
    heat.promoter.data$baseMean
  ),
  method = "spearman",
  exact = FALSE
)
cor.mean.expression
# Spearman's rank correlation rho
# 
# data:  log10(heat.promoter.data$mean_normalized_CPM) and log10(heat.promoter.data$baseMean)
# S = 1745815, p-value = 1.988e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2699747 

# 59. Promoter accessibility vs magnitude of heat response


cor.heat.response <- cor.test(
  log10(
    heat.promoter.data$mean_normalized_CPM
  ),
  abs(
    heat.promoter.data$log2FoldChange
  ),
  method = "spearman",
  exact = FALSE
)

cor.heat.response

# Spearman's rank correlation rho
# 
# data:  log10(heat.promoter.data$mean_normalized_CPM) and abs(heat.promoter.data$log2FoldChange)
# S = 2408714, p-value = 0.9108
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#          rho 
# -0.007221583


# 60. Save correlation statistics

correlation.summary <- data.frame(
  comparison = c(
    "Promoter accessibility vs mean expression",
    "Promoter accessibility vs absolute heat log2FC"
  ),
  n_genes = c(
    nrow(heat.promoter.data),
    nrow(heat.promoter.data)
  ),
  spearman_rho = c(
    unname(
      cor.mean.expression$estimate
    ),
    unname(
      cor.heat.response$estimate
    )
  ),
  p_value = c(
    cor.mean.expression$p.value,
    cor.heat.response$p.value
  )
)


write.csv(
  correlation.summary,
  file = file.path(
    output.dir,
    "ATAC_RNA_Spearman_correlations.csv"
  ),
  row.names = FALSE
)



# PART M: Figure 4 statistical labels

# 61. Format mean-expression p-value label

if (cor.mean.expression$p.value < 0.001) {
  
  mean.expression.p.label <- "p < 0.001"
  
} else {
  
  mean.expression.p.label <- paste0(
    "p = ",
    formatC(
      cor.mean.expression$p.value,
      format = "f",
      digits = 3
    )
  )
  
}


mean.expression.label <- paste0(
  "Spearman rho = ",
  sprintf(
    "%.2f",
    unname(
      cor.mean.expression$estimate
    )
  ),
  "\n",
  mean.expression.p.label
)


# 63. Format heat-response p-value label

if (cor.heat.response$p.value < 0.001) {
  
  heat.response.p.label <- "p < 0.001"
  
} else {
  
  heat.response.p.label <- paste0(
    "p = ",
    formatC(
      cor.heat.response$p.value,
      format = "f",
      digits = 3
    )
  )
  
}


heat.response.label <- paste0(
  "Spearman rho = ",
  sprintf(
    "%.2f",
    unname(
      cor.heat.response$estimate
    )
  ),
  "\n",
  heat.response.p.label
)


# PART N: Figure 4

# 64. Define shared Figure 4 theme


figure4.theme <- theme_pubr(
  base_size = 22
) +
  theme(
    plot.title = element_text(
      size = 24,
      face = "plain",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 22,
      face = "plain"
    ),
    axis.text = element_text(
      size = 18
    ),
    plot.margin = margin(
      15,
      20,
      15,
      15
    )
  )



# 64. Figure 4A:
#     promoter accessibility vs mean expression

plot.access.mean.expression <- ggplot(
  heat.promoter.data,
  aes(
    x = log10(mean_normalized_CPM),
    y = log10(baseMean)
  )
) +
  geom_point(
    alpha = 0.6,
    size = 3,
    color = "grey35"
  ) +
  geom_smooth(
    method = "lm",
    color = "black",
    se = TRUE,
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = mean.expression.label,
    hjust = 1.05,
    vjust = 1.15,
    size = 7,
    fontface = "italic"
  ) +
  labs(
    x = expression(
      log[10] * " promoter accessibility (mean normalized CPM)"
    ),
    y = expression(
      log[10] * " mean normalized expression"
    )
  ) +
  figure4.theme


plot.access.mean.expression

# 66. Figure 4B:
#     promoter accessibility vs magnitude of heat response


plot.access.heat.response <- ggplot(
  heat.promoter.data,
  aes(
    x = log10(mean_normalized_CPM),
    y = abs(log2FoldChange)
  )
) +
  geom_point(
    alpha = 0.6,
    size = 3,
    color = "grey35"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = heat.response.label,
    hjust = 1.05,
    vjust = 1.15,
    size = 7,
    fontface = "italic"
  ) +
  labs(
    x = expression(
      log[10] * " promoter accessibility (mean normalized CPM)"
    ),
    y = expression(
      "abs. value log"[2] * " fold change"
    ),
  ) +
  figure4.theme


plot.access.heat.response



# 67. Combine Figure 4 panels


figure4 <- plot.access.mean.expression +
  plot.access.heat.response +
  plot_annotation(
    tag_levels = "A"
  )



figure4



# 68. Save Figure 4


ggsave(
  filename = file.path(
    output.dir,
    "Figure4_ATAC_heat_response_correlations.pdf"
  ),
  plot = figure4,
  width = 16,
  height = 7
)


ggsave(
  filename = file.path(
    output.dir,
    "Figure4_ATAC_heat_response_correlations.tiff"
  ),
  plot = figure4,
  width = 14,
  height = 7,
  dpi = 300,
  compression = "lzw"
)



# PART N: Save Figure 3 analysis tables


# 69. Save expression-quintile data


write.csv(
  promoter.expression.data,
  file = file.path(
    output.dir,
    "Figure3A_expression_quintile_data.csv"
  ),
  row.names = FALSE
)



# 70. Save expression-variability data

write.csv(
  promoter.variation.extremes,
  file = file.path(
    output.dir,
    "Figure3B_expression_variability_data.csv"
  ),
  row.names = FALSE
)



# PART O: Final analysis summary

# 71. Print final summary




cat(
  "RNA-seq expression universe:",
  length(expression.universe),
  "\n"
)


cat(
  "Expressed genes with accessible promoters:",
  nrow(closest.promoter.per.gene),
  "\n\n"
)


cat(
  "Heat-responsive genes:",
  length(heat.genes),
  "\n"
)


cat(
  "Heat-responsive genes with accessible promoters:",
  heat.accessible,
  "\n"
)


cat(
  "Heat-response Fisher odds ratio:",
  round(
    unname(
      heat.promoter.fisher$estimate
    ),
    3
  ),
  "\n"
)


cat(
  "Heat-response Fisher p-value:",
  signif(
    heat.promoter.fisher$p.value,
    4
  ),
  "\n\n"
)


cat(
  "Dixon red-module genes represented:",
  length(dixon.red.universe),
  "\n"
)


cat(
  "Dixon red-module genes with accessible promoters:",
  dixon.accessible,
  "\n"
)


cat(
  "Dixon Fisher odds ratio:",
  round(
    unname(
      dixon.promoter.fisher$estimate
    ),
    3
  ),
  "\n"
)


cat(
  "Dixon Fisher p-value:",
  signif(
    dixon.promoter.fisher$p.value,
    4
  ),
  "\n\n"
)


cat(
  "Heat-responsive genes used for Figure 4:",
  nrow(heat.promoter.data),
  "\n"
)


cat(
  "Accessibility vs mean expression rho:",
  round(
    unname(
      cor.mean.expression$estimate
    ),
    3
  ),
  "\n"
)


cat(
  "Accessibility vs mean expression p:",
  signif(
    cor.mean.expression$p.value,
    4
  ),
  "\n"
)


cat(
  "Accessibility vs |heat log2FC| rho:",
  round(
    unname(
      cor.heat.response$estimate
    ),
    3
  ),
  "\n"
)


cat(
  "Accessibility vs |heat log2FC| p:",
  signif(
    cor.heat.response$p.value,
    4
  ),
  "\n"
)



# 72. Save plain-text summary
capture.output(
  {
    
    cat(
      "ATAC-seq / RNA-seq integration\n"
    )
    
    cat(
      "==============================\n\n"
    )
    
    cat(
      "Interpretation:\n"
    )
    
    cat(
      "ATAC-seq and RNA-seq datasets were generated from\n"
    )
    
    cat(
      "independently sourced fragments. Results therefore\n"
    )
    
    cat(
      "represent broad cross-dataset associations.\n\n"
    )
    
    cat(
      "RNA-seq expression universe:",
      length(expression.universe),
      "\n"
    )
    
    cat(
      "Expressed genes with accessible promoters:",
      nrow(closest.promoter.per.gene),
      "\n\n"
    )
    
    cat(
      "Heat-responsive promoter enrichment\n"
    )
    
    cat(
      "-----------------------------------\n"
    )
    
    cat(
      "Heat-responsive genes:",
      length(heat.genes),
      "\n"
    )
    
    cat(
      "Heat-responsive genes with accessible promoters:",
      heat.accessible,
      "\n"
    )
    
    cat(
      "Odds ratio:",
      unname(
        heat.promoter.fisher$estimate
      ),
      "\n"
    )
    
    cat(
      "95% CI:",
      heat.promoter.fisher$conf.int[1],
      "to",
      heat.promoter.fisher$conf.int[2],
      "\n"
    )
    
    cat(
      "p-value:",
      heat.promoter.fisher$p.value,
      "\n\n"
    )
    
    cat(
      "Dixon red-module promoter enrichment\n"
    )
    
    cat(
      "------------------------------------\n"
    )
    
    cat(
      "Red-module genes represented:",
      length(dixon.red.universe),
      "\n"
    )
    
    cat(
      "Red-module genes with accessible promoters:",
      dixon.accessible,
      "\n"
    )
    
    cat(
      "Odds ratio:",
      unname(
        dixon.promoter.fisher$estimate
      ),
      "\n"
    )
    
    cat(
      "95% CI:",
      dixon.promoter.fisher$conf.int[1],
      "to",
      dixon.promoter.fisher$conf.int[2],
      "\n"
    )
    
    cat(
      "p-value:",
      dixon.promoter.fisher$p.value,
      "\n\n"
    )
    
    cat(
      "Figure 4 correlations\n"
    )
    
    cat(
      "---------------------\n"
    )
    
    cat(
      "Genes used:",
      nrow(heat.promoter.data),
      "\n"
    )
    
    cat(
      "Accessibility vs mean expression:\n"
    )
    
    cat(
      "Spearman rho:",
      unname(
        cor.mean.expression$estimate
      ),
      "\n"
    )
    
    cat(
      "p-value:",
      cor.mean.expression$p.value,
      "\n\n"
    )
    
    cat(
      "Accessibility vs magnitude of heat response:\n"
    )
    
    cat(
      "Spearman rho:",
      unname(
        cor.heat.response$estimate
      ),
      "\n"
    )
    
    cat(
      "p-value:",
      cor.heat.response$p.value,
      "\n"
    )
    
  },
  file = file.path(
    output.dir,
    "ATAC_RNA_integration_summary.txt"
  )
)



# 73. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "ATAC_RNA_integration_sessionInfo.txt"
  )
)

