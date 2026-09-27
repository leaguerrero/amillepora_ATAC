# 07_ATAC_QC.R
#
# Purpose:
#   Summarize quality-control characteristics of the four
#   Acropora millepora ATAC-seq libraries using the counted
#   DiffBind object generated in Script 06.
#
#   This script:
#     1. Loads the counted DiffBind object
#     2. Confirms the four ATAC-seq samples
#     3. Extracts DiffBind sample information
#     4. Reports DiffBind read counts and FRiP
#     5. Calculates estimated reads within peaks
#     6. Generates the DiffBind correlation plot
#     7. Generates the DiffBind PCA plot
#
# Inputs from Script 06:
#   - ATAC_DiffBind_object.rds
#
# Outputs:
#   - ATAC_DiffBind_sample_information.csv
#   - ATAC_DiffBind_FRiP.csv
#   - ATAC_DiffBind_correlation_plot.pdf
#   - ATAC_DiffBind_correlation_plot.tiff
#   - ATAC_DiffBind_PCA.pdf
#   - ATAC_DiffBind_PCA.tiff
#   - ATAC_QC_summary.txt
#   - ATAC_QC_sessionInfo.txt
#
# Important:
#   This script summarizes QC information contained in the
#   DiffBind object. 



# 1. Load package

library(DiffBind)



# 2. Define input and output directories

project.dir <- "./Code_and_data/"

# Counted DiffBind object generated in Script 06

diffbind.file <- file.path(
  project.dir,
  "./Data/06_consensus_peaks",
  "ATAC_DiffBind_object.rds"
)


# Output directory

output.dir <- file.path(
  project.dir,
  "./Data/07_QC"
)


dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# PART A: Load and inspect DiffBind object

# 3. Load counted DiffBind object


DBdata <- readRDS(
  diffbind.file
)


cat(
  "Counted DiffBind object loaded.\n"
)



# 4. Display DiffBind object


cat(
  "\nDiffBind object summary:\n"
)


print(
  DBdata
)

# 4 Samples, 11231 sites in matrix:
#   ID    Reads FRiP
# 1  1  3066830 0.11
# 2  2 15199132 0.03
# 3  3  5544314 0.02
# 4  4 12873484 0.02

# 5. Extract sample information from DiffBind


sample.info <- dba.show(
  DBdata
)


cat(
  "\nDiffBind sample information:\n"
)


print(
  sample.info
)

# ID Caller Intervals    Reads FRiP
# 1  1 counts     11231  3066830 0.11
# 2  2 counts     11231 15199132 0.03
# 3  3 counts     11231  5544314 0.02
# 4  4 counts     11231 12873484 0.02

cat(
  "\nNumber of ATAC-seq samples:",
  nrow(sample.info),
  "\n"
)


# Expected sample number

if (nrow(sample.info) != 4) {
  
  warning(
    "Expected 4 ATAC-seq samples."
  )
  
}

# 6. Save complete DiffBind sample information
write.csv(
  sample.info,
  file = file.path(
    output.dir,
    "ATAC_DiffBind_sample_information.csv"
  ),
  row.names = FALSE
)

# PART B: Read depth and FRiP
# 7. Confirm that required QC columns are present


if (!"Reads" %in% colnames(sample.info)) {
  
  stop(
    "The DiffBind sample information does not contain a Reads column."
  )
  
}


if (!"FRiP" %in% colnames(sample.info)) {
  
  stop(
    "The DiffBind sample information does not contain a FRiP column."
  )
  
}



# 8. Extract sample IDs
# The original analysis used the DiffBind ID column as the
# sample identifier.


if ("ID" %in% colnames(sample.info)) {
  
  sample.id <- sample.info$ID
  
} else {
  
  sample.id <- rownames(sample.info)
  
}


cat(
  "\nSample IDs:\n"
)


print(
  sample.id
)

# [1] "1" "2" "3" "4"

# 9. Extract DiffBind read counts

# These values come directly from dba.show().
#
# They should not be confused with the total number of raw
# sequencing reads generated before filtering and processing.


diffbind.reads <- sample.info$Reads


cat(
  "\nDiffBind Reads:\n"
)


print(
  diffbind.reads
)

# [1]  3066830 15199132  5544314 12873484

# 10. Extract FRiP values
# FRiP = Fraction of Reads in Peaks


frip <- sample.info$FRiP


cat(
  "\nFRiP values:\n"
)


print(
  frip
)

# [1] 0.11 0.03 0.02 0.02

# 11. Estimate number of reads within peaks


# This reproduces the calculation in the original analysis:
#     PeakReads = Reads x FRiP


peak.reads <- round(
  diffbind.reads * frip
)


cat(
  "\nEstimated reads within peaks:\n"
)


print(
  peak.reads
)

# [1] 337351 455974 110886 257470

# 12. Build QC table


frip.table <- data.frame(
  sample = sample.id,
  DiffBind_Reads = diffbind.reads,
  FRiP = frip,
  Estimated_Peak_Reads = peak.reads,
  stringsAsFactors = FALSE
)


cat(
  "\nATAC-seq DiffBind QC table:\n"
)


print(
  frip.table
)
# sample DiffBind_Reads FRiP Estimated_Peak_Reads
# 1      1        3066830 0.11               337351
# 2      2       15199132 0.03               455974
# 3      3        5544314 0.02               110886
# 4      4       12873484 0.02               257470


# 13. Save QC table


write.csv(
  frip.table,
  file = file.path(
    output.dir,
    "ATAC_DiffBind_FRiP.csv"
  ),
  row.names = FALSE
)



# 14. Summarize FRiP values


mean.frip <- mean(
  frip,
  na.rm = TRUE
)


minimum.frip <- min(
  frip,
  na.rm = TRUE
)


maximum.frip <- max(
  frip,
  na.rm = TRUE
)


cat(
  "\nMean FRiP:",
  round(mean.frip, 4),
  "\n"
)

# Mean FRiP: 0.045 

cat(
  "Minimum FRiP:",
  round(minimum.frip, 4),
  "\n"
)

# Minimum FRiP: 0.02 

cat(
  "Maximum FRiP:",
  round(maximum.frip, 4),
  "\n"
)
# Maximum FRiP: 0.11 


# 15. Summarize DiffBind read counts

mean.diffbind.reads <- mean(
  diffbind.reads,
  na.rm = TRUE
)


minimum.diffbind.reads <- min(
  diffbind.reads,
  na.rm = TRUE
)


maximum.diffbind.reads <- max(
  diffbind.reads,
  na.rm = TRUE
)


cat(
  "\nMean DiffBind Reads:",
  round(mean.diffbind.reads),
  "\n"
)
# Mean DiffBind Reads: 9170940 

cat(
  "Minimum DiffBind Reads:",
  minimum.diffbind.reads,
  "\n"
)
# Minimum DiffBind Reads: 3066830 


cat(
  "Maximum DiffBind Reads:",
  maximum.diffbind.reads,
  "\n"
)
# Maximum DiffBind Reads: 15199132 


# PART C: DiffBind sample correlation

# 16. Display DiffBind correlation plot


plot(
  DBdata
)



# 17. Save DiffBind correlation plot as PDF


pdf(
  file = file.path(
    output.dir,
    "ATAC_DiffBind_correlation_plot.pdf"
  ),
  width = 8,
  height = 8
)


plot(
  DBdata
)


dev.off()



# 18. Save DiffBind correlation plot as TIFF


tiff(
  filename = file.path(
    output.dir,
    "ATAC_DiffBind_correlation_plot.tiff"
  ),
  width = 8,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)


plot(
  DBdata
)


dev.off()



# PART D: ATAC-seq PCA




# 19. Display DiffBind PCA
# This reproduces the PCA call from the original analysis:
#
#     dba.plotPCA(DBdata, label = DBA_ID)
#
# Each point represents one of the four independent
# ATAC-seq samples.

dba.plotPCA(
  DBdata,
  label = DBA_ID
)



# 20. Save ATAC-seq PCA as PDF


pdf(
  file = file.path(
    output.dir,
    "ATAC_DiffBind_PCA.pdf"
  ),
  width = 8,
  height = 7
)


dba.plotPCA(
  DBdata,
  label = DBA_ID
)


dev.off()



# 21. Save ATAC-seq PCA as TIFF


tiff(
  filename = file.path(
    output.dir,
    "ATAC_DiffBind_PCA.tiff"
  ),
  width = 8,
  height = 7,
  units = "in",
  res = 300,
  compression = "lzw"
)


dba.plotPCA(
  DBdata,
  label = DBA_ID
)


dev.off()



# PART E: Final QC summary

# 22. Print final summary



cat(
  "Samples:",
  nrow(sample.info),
  "\n"
)


cat(
  "\nDiffBind read counts:\n"
)


print(
  diffbind.reads
)


cat(
  "\nFRiP:\n"
)


print(
  frip
)


cat(
  "\nMean FRiP:",
  round(mean.frip, 4),
  "\n"
)


cat(
  "FRiP range:",
  round(minimum.frip, 4),
  "to",
  round(maximum.frip, 4),
  "\n"
)



# 23. Save plain-text QC summary


capture.output(
  {
    
    cat(
      "ATAC-seq DiffBind quality control\n"
    )
    
    
    cat(
      "Number of ATAC-seq samples:",
      nrow(sample.info),
      "\n\n"
    )
    
    cat(
      "DiffBind sample information:\n"
    )
    
    print(
      sample.info
    )
    
    cat(
      "\nFRiP summary:\n"
    )
    
    print(
      frip.table
    )
    
    cat(
      "\nMean FRiP:",
      round(mean.frip, 4),
      "\n"
    )
    
    cat(
      "Minimum FRiP:",
      round(minimum.frip, 4),
      "\n"
    )
    
    cat(
      "Maximum FRiP:",
      round(maximum.frip, 4),
      "\n\n"
    )
    
    cat(
      "Mean DiffBind Reads:",
      round(mean.diffbind.reads),
      "\n"
    )
    
    cat(
      "Minimum DiffBind Reads:",
      minimum.diffbind.reads,
      "\n"
    )
    
    cat(
      "Maximum DiffBind Reads:",
      maximum.diffbind.reads,
      "\n\n"
    )
    
    cat(
      "Additional ATAC-seq QC metrics described in the manuscript\n"
    )
    
    cat(
      "were generated during upstream processing and are not\n"
    )
    
    cat(
      "recalculated in this R script:\n"
    )
    
    cat(
      "  - Bioanalyzer fragment-size profiles\n"
    )
    
    cat(
      "  - library complexity\n"
    )

    
    cat(
      "  - TSS enrichment\n"
    )
    
  },
  file = file.path(
    output.dir,
    "ATAC_QC_summary.txt"
  )
)



# 24. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "ATAC_QC_sessionInfo.txt"
  )
)
