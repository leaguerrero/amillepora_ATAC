# 01_RNAseq_build_count_matrix.R
#
# Purpose:
#   Read HTSeq gene-count files and sample metadata.
#   Construct the RNA-seq count matrix used for downstream analyses.
#   Remove genes with fewer than 10 total reads across all samples.
#
# Inputs:
#   - HTSeq count files
#   - tag_seq_meta_data.txt
#
# Outputs:
#   - RNAseq_sample_table.csv
#   - RNAseq_sample_to_HTSeq_file_mapping.csv
#   - RNAseq_raw_count_matrix.csv
#   - RNAseq_filtered_count_matrix.csv
#
# Expected checkpoints:
#   - 58 RNA-seq samples
#   - 19,649 genes retained after filtering


# 1. Load package
library(DESeq2)


# 2. Define input and output directories


project.dir <- "./Code_and_data/"

count.dir <- file.path(
  project.dir,
  "./Data/00_inputs/RNAseq_counts/"
)

metadata.file <- file.path(
  project.dir,
  "./Data/00_inputs/RNAseq_counts/tag_seq_meta_data.txt"
)

output.dir <- file.path(
  project.dir,
  "./Data/01_count_matrix"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Identify HTSeq count files

count.files <- list.files(
  count.dir,
  pattern = "\\gene_count.txt$",
  full.names = FALSE
)

count.files <- sort(count.files)

cat(
  "Number of HTSeq count files:",
  length(count.files),
  "\n"
)

# Number of HTSeq count files: 58 


# 4. Read sample metadata

metadata <- read.table(
  metadata.file,
  header = TRUE,
  stringsAsFactors = FALSE
)

cat(
  "Number of rows in metadata:",
  nrow(metadata),
  "\n"
)

# Number of rows in metadata: 58 


print(
  colnames(metadata)
)
# [1] "ID"                "acclimation.temp"  "acclimation.class"
# [4] "acclimation.txt"   "stress.txt"        "stress.group"     
# [7] "stress.day"        "ext.date"          "ext.batch"        
# [10] "genotype"          "file.name"  


# 5. Define sample names

sample.names <- c(
  "01_27_Y", "01_28_R", "01_40_Y", "01_42_R", "01_45_Y",
  "01_46_Y", "01_51_Y", "01_52_Y", "01_54_Y", "01_58_Y",
  "01_63_Y", "01_68_Y", "01_78_R", "01_92_Y", "02_29_Y",
  "02_34_Y", "02_35_R", "02_47_Y", "02_65_R", "02_79_R",
  "03_29_R", "03_48_Y", "03_62_R", "03_72_R", "03_90_Y",
  "03_91_R", "04_48_R", "04_65_Y", "04_70_Y", "04_95_Y",
  "05_43_R", "05_46_R", "05_50_Y", "05_62_Y", "05_82_Y",
  "05_88_R", "06_22_Y", "06_32_R", "06_59_R", "06_64_R",
  "06_67_R", "06_72_Y", "06_82_R", "06_92_R", "07_31_R",
  "07_36_R", "07_39_R", "07_53_Y", "07_64_Y", "07_88_Y",
  "08_19_R", "08_30_R", "08_33_Y", "08_36_Y", "08_45_R",
  "08_50_R", "08_85_Y", "08_95_R"
)

cat(
  "Number of sample names:",
  length(sample.names),
  "\n"
)

# Number of sample names: 58 



# 6. Check correspondence between samples and HTSeq files

# Inspect table to confirm that each sample name is associated
# with the correct HTSeq count file.


sample.file.mapping <- data.frame(
  sampleName = sample.names,
  fileName = count.files
)

print(
  sample.file.mapping,
  row.names = FALSE
)

write.csv(
  sample.file.mapping,
  file = file.path(
    output.dir,
    "RNAseq_sample_to_HTSeq_file_mapping.csv"
  ),
  row.names = FALSE
)



# 7. Construct sample table for DESeq2

# Columns 2 through 10 are the experimental metadata columns
# used in the original analysis.


sample.table <- data.frame(
  sampleName = sample.names,
  fileName = count.files,
  metadata[, 2:10],
  stringsAsFactors = FALSE
)

rownames(sample.table) <- sample.table$sampleName

cat(
  "\nRNA-seq sample table:\n"
)

print(
  sample.table
)


# Save the sample table used to construct the count matrix

write.csv(
  sample.table,
  file = file.path(
    output.dir,
    "RNAseq_sample_table.csv"
  ),
  row.names = FALSE
)



# 8. Construct DESeq2 dataset from HTSeq files

# No statistical model is being fit here.
# The design is therefore set to ~ 1.


dds.counts <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample.table,
  directory = count.dir,
  design = ~ 1
)



# 9. Extract raw gene-count matrix


raw.counts <- counts(
  dds.counts,
  normalized = FALSE
)

cat(
  "\nNumber of genes before filtering:",
  nrow(raw.counts),
  "\n"
)

# Number of genes before filtering: 38247 

cat(
  "Number of samples:",
  ncol(raw.counts),
  "\n"
)

# Number of samples: 58 

# Add gene IDs as the first column for the output file

raw.counts.output <- data.frame(
  gene = rownames(raw.counts),
  raw.counts,
  check.names = FALSE
)

write.csv(
  raw.counts.output,
  file = file.path(
    output.dir,
    "RNAseq_raw_count_matrix.csv"
  ),
  row.names = FALSE
)



# 11. Filter low-count genes

# Retain genes with at least 10 total reads across all 58 samples.


keep.genes <- rowSums(raw.counts) >= 10

filtered.counts <- raw.counts[
  keep.genes,
  ,
  drop = FALSE
]

cat(
  "\nGenes retained after filtering:",
  nrow(filtered.counts),
  "\n"
)

# Genes retained after filtering: 19649 

cat(
  "Genes removed by filtering:",
  sum(!keep.genes),
  "\n"
)

# Genes removed by filtering: 18598 

# 12. Save filtered count matrix


filtered.counts.output <- data.frame(
  gene = rownames(filtered.counts),
  filtered.counts,
  check.names = FALSE
)

write.csv(
  filtered.counts.output,
  file = file.path(
    output.dir,
    "RNAseq_filtered_count_matrix.csv"
  ),
  row.names = FALSE
)



# 13. Final checks


cat(
  "\n========================================\n"
)

cat(
  "RNA-seq count-matrix construction complete\n"
)

cat(
  "========================================\n"
)

cat(
  "Samples:",
  ncol(filtered.counts),
  "\n"
)

cat(
  "Genes before filtering:",
  nrow(raw.counts),
  "\n"
)

cat(
  "Genes after filtering:",
  nrow(filtered.counts),
  "\n"
)


# Expected number of samples

if (ncol(filtered.counts) != 58) {
  warning(
    "Expected 58 samples. Check sample and HTSeq file mapping."
  )
}


# Expected number of genes based on the final manuscript analysis

if (nrow(filtered.counts) != 19649) {
  warning(
    "Expected 19,649 genes after filtering. Check the input files and filtering."
  )
}



# 14. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "RNAseq_count_matrix_sessionInfo.txt"
  )
)

