# 06_ATAC_build_consensus_peaks.R
#
# Purpose:
#   Construct the consensus ATAC-seq peak set used for
#   downstream analyses in Acropora millepora.
#
#   This script:
#     1. Reads the DiffBind sample sheet
#     2. Creates a DiffBind object
#     3. Counts reads across binding intervals
#     4. Identifies consensus peaks present in at least
#        2 of the 4 ATAC-seq samples
#     5. Applies the MACS peak-confidence threshold
#     6. Normalizes ATAC-seq read counts using CPM
#     7. Calculates mean accessibility across the four samples
#     8. Annotates consensus peaks using ChIPseeker
#     9. Saves processed ATAC-seq data for downstream analyses
#
# Inputs:
#   - sampleSheet.csv
#   - peak files referenced by sampleSheet.csv
#   - BAM files referenced by sampleSheet.csv
#   - Amil.all.maker.noseq.gff
#
# Outputs:
#   - ATAC_DiffBind_object.rds
#   - ATAC_consensus_peaks.rds
#   - ATAC_consensus_peaks.csv
#   - ATAC_consensus_peaks_2of4.bed
#   - ATAC_consensus_peak_counts.csv
#   - ATAC_consensus_peak_CPM.csv
#   - ATAC_annotated_consensus_peaks.csv
#   - ATAC_unannotated_consensus_peaks.csv
#   - ATAC_consensus_peak_summary.txt
#
# Expected checkpoints:
#   - 4 ATAC-seq samples
#   - 11,231 consensus peaks
#   - 3,973 promoter-associated peaks




# 1. Load packages
library(DiffBind)
library(edgeR)
library(GenomicFeatures)
library(ChIPseeker)



# 2. Define input and output paths


project.dir <- "./Code_and_data/"


# DiffBind sample sheet

sample.sheet.file <- file.path(
  project.dir,
  "Data/00_inputs/ATAC_peaks/sampleSheet.csv"
)



# A. millepora gene annotation

gff.file <- file.path(
  project.dir,
  "./external_resources/Amil_v2.01",
  "Amil.all.maker.noseq.gff"
)


# Output directory

output.dir <- file.path(
  project.dir,
  "./Data/06_consensus_peaks"
)


dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# PART A: Build DiffBind object
# 3. Read ATAC-seq sample sheet

samples <- read.csv(
  sample.sheet.file,
  header = TRUE,
  stringsAsFactors = FALSE
)


cat(
  "Number of samples in ATAC-seq sample sheet:",
  nrow(samples),
  "\n"
)

# Number of samples in ATAC-seq sample sheet: 4 

cat(
  "\nATAC-seq sample sheet:\n"
)

print(
  samples
)


# SampleID
# 1        1
# 2        2
# 3        3
# 4        4
# bamReads
# 1 ./Code_and_data/Data/00_inputs/ATAC_peaks/BAM_files/AMIL_GENO_01_S1_nuc_deduplicated_reads.rg.bam
# 2 ./Code_and_data/Data/00_inputs/ATAC_peaks/BAM_files/AMIL_GENO_02_S2_nuc_deduplicated_reads.rg.bam
# 3 ./Code_and_data/Data/00_inputs/ATAC_peaks/BAM_files/AMIL_GENO_03_S3_nuc_deduplicated_reads.rg.bam
# 4 ./Code_and_data/Data/00_inputs/ATAC_peaks/BAM_files/AMIL_GENO_04_S4_nuc_deduplicated_reads.rg.bam
# Peaks
# 1 ./Code_and_data/Data/00_inputs/ATAC_peaks/Peak_files/AMIL_GENO_01_S1_peaks.xls
# 2 ./Code_and_data/Data/00_inputs/ATAC_peaks/Peak_files/AMIL_GENO_02_S2_peaks.xls
# 3 ./Code_and_data/Data/00_inputs/ATAC_peaks/Peak_files/AMIL_GENO_03_S3_peaks.xls
# 4 ./Code_and_data/Data/00_inputs/ATAC_peaks/Peak_files/AMIL_GENO_04_S4_peaks.xls
# PeakCaller
# 1       macs
# 2       macs
# 3       macs
# 4       macs

# Expected number of samples

if (nrow(samples) != 4) {
  
  warning(
    "Expected 4 ATAC-seq samples."
  )
  
}

# 4. Create DiffBind object


DBdata <- dba(
  sampleSheet = samples
)


cat(
  "\nInitial DiffBind object:\n"
)

print(
  DBdata
)

# 4 Samples, 19400 sites in matrix (97847 total):
#   ID Intervals
# 1  1     11777
# 2  2     81283
# 3  3     10708
# 4  4     25055

# 5. Count reads across binding-site intervals


DBdata <- dba.count(
  DBdata
)


cat(
  "\nDiffBind object after counting reads:\n"
)

print(
  DBdata
)



# 6. Save counted DiffBind object
# This preserves the counted ATAC-seq dataset so later scripts
# do not need to repeat read counting.

saveRDS(
  DBdata,
  file = file.path(
    output.dir,
    "ATAC_DiffBind_object.rds"
  )
)

# PART B: Construct consensus peak set



# 7. Extract consensus peaks
# A consensus peak must:
#
#   1. Be detected in at least 2 of the 4 samples
#
#   2. Pass a MACS -log10(p-value) score threshold of 4.3
#
# The 4.3 threshold corresponds to:
#
#   p = 5 x 10^-5



consensus.peaks <- dba.peakset(
  DBdata,
  minOverlap = 2,
  peak.caller = "macs",
  peak.format = "macs",
  filter = 4.3,
  scoreCol = 7,
  bRetrieve = TRUE
)


cat(
  "\nTotal consensus peaks:",
  length(consensus.peaks),
  "\n"
)

# Total consensus peaks: 11231 

# Expected checkpoint

if (length(consensus.peaks) != 11231) {
  
  warning(
    "Expected 11,231 consensus ATAC-seq peaks."
  )
  
}

# 8. Create unique coordinate identifier for each peak

consensus.peak.id <- paste0(
  seqnames(consensus.peaks),
  ":",
  start(consensus.peaks),
  "-",
  end(consensus.peaks)
)


# Confirm that every consensus peak has a unique coordinate

if (anyDuplicated(consensus.peak.id) > 0) {
  
  stop(
    "Duplicate consensus peak coordinates were detected."
  )
  
}


cat(
  "Unique consensus peak coordinates:",
  length(unique(consensus.peak.id)),
  "\n"
)
# Unique consensus peak coordinates: 11231 


# 9. Save consensus peak GRanges object


saveRDS(
  consensus.peaks,
  file = file.path(
    output.dir,
    "ATAC_consensus_peaks.rds"
  )
)

# 10. Save consensus peak coordinates as CSV


consensus.peaks.table <- as.data.frame(
  consensus.peaks
)


consensus.peaks.table$peak_id <- consensus.peak.id


# Put peak identifier first

consensus.peaks.table <- consensus.peaks.table[
  ,
  c(
    "peak_id",
    setdiff(
      colnames(consensus.peaks.table),
      "peak_id"
    )
  )
]


write.csv(
  consensus.peaks.table,
  file = file.path(
    output.dir,
    "ATAC_consensus_peaks.csv"
  ),
  row.names = FALSE
)



# 11. Save consensus peaks in BED format
# GRanges uses 1-based coordinates.
# BED format uses a 0-based start coordinate.
#
# Therefore:
#
#   BED start = GRanges start - 1
#   BED end   = GRanges end

consensus.peaks.bed <- data.frame(
  chromosome = as.character(
    seqnames(consensus.peaks)
  ),
  start = start(consensus.peaks) - 1,
  end = end(consensus.peaks),
  peak_id = consensus.peak.id,
  stringsAsFactors = FALSE
)


write.table(
  consensus.peaks.bed,
  file = file.path(
    output.dir,
    "ATAC_consensus_peaks_2of4.bed"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



# PART C: Normalize ATAC-seq accessibility

# 12. Extract DiffBind read-count matrix
# This is the same count matrix used in the original analysis.

atac.counts <- DBdata$binding


cat(
  "\nATAC-seq count-matrix dimensions:\n"
)

print(
  dim(atac.counts)
)
# [1] 11231     7

cat(
  "\nATAC-seq count-matrix column names:\n"
)

print(
  colnames(atac.counts)
)
# 1       2       3       4 
# "CHR" "START"   "END"     "1"     "2"     "3"     "4" 


# 13. Confirm one count-matrix row per consensus peak
# The original analysis attached row-wise CPM values from
# DBdata$binding directly to the consensus peak set.


if (nrow(atac.counts) != length(consensus.peaks)) {
  
  stop(
    paste(
      "ATAC count matrix contains",
      nrow(atac.counts),
      "rows but the consensus peak set contains",
      length(consensus.peaks),
      "peaks.",
      "Stop here and inspect the DiffBind object before continuing."
    )
  )
  
}


cat(
  "ATAC count matrix contains one row per consensus peak.\n"
)

# ATAC count matrix contains one row per consensus peak.

# 14. Normalize ATAC-seq counts using CPM


atac.cpm <- cpm(
  atac.counts
)


cat(
  "\nCPM-normalized matrix dimensions:\n"
)

print(
  dim(atac.cpm)
)
# [1] 11231     7


# 15. Calculate mean accessibility across the four samples


avg.accessibility <- rowMeans(
  atac.cpm
)


cat(
  "\nSummary of mean CPM accessibility:\n"
)

print(
  summary(avg.accessibility)
)

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 3.45    36.76    61.89    89.04    96.61 92628.99 


# 16. Add mean accessibility to consensus peaks


mcols(consensus.peaks)$avgAccessibility <- avg.accessibility



# 17. Save raw ATAC-seq count matrix

atac.counts.output <- data.frame(
  peak_id = consensus.peak.id,
  atac.counts,
  check.names = FALSE
)


write.csv(
  atac.counts.output,
  file = file.path(
    output.dir,
    "ATAC_consensus_peak_counts.csv"
  ),
  row.names = FALSE
)



# 18. Save CPM-normalized ATAC-seq matrix

atac.cpm.output <- data.frame(
  peak_id = consensus.peak.id,
  atac.cpm,
  mean_CPM = avg.accessibility,
  check.names = FALSE
)


write.csv(
  atac.cpm.output,
  file = file.path(
    output.dir,
    "ATAC_consensus_peak_CPM.csv"
  ),
  row.names = FALSE
)



# PART D: Annotate consensus peaks

# 19. Create transcript database from A. millepora GFF


txdb <- makeTxDbFromGFF(
  gff.file,
  format = "gff3"
)



# 20. Annotate consensus peaks with ChIPseeker


consensus.peaks.anno <- annotatePeak(
  consensus.peaks,
  TxDb = txdb
)

# 21. Convert ChIPseeker annotation to data frame


consensus.peaks.anno.df <- as.data.frame(
  consensus.peaks.anno
)


# 22. Create coordinate IDs for annotated peaks


annotated.peak.id <- paste0(
  consensus.peaks.anno.df$seqnames,
  ":",
  consensus.peaks.anno.df$start,
  "-",
  consensus.peaks.anno.df$end
)


consensus.peaks.anno.df$peak_id <- annotated.peak.id



# 23. Match accessibility values to annotated peaks

# ChIPseeker may return fewer peaks than were supplied.
#
# Accessibility is therefore matched back to the consensus
# peak set by genomic coordinate rather than assuming that
# row numbers are identical.


match.index <- match(
  annotated.peak.id,
  consensus.peak.id
)


# Confirm every annotated peak matches a consensus peak

if (any(is.na(match.index))) {
  
  stop(
    "At least one annotated peak could not be matched to the consensus peak set."
  )
  
}


consensus.peaks.anno.df$avgAccess <- avg.accessibility[
  match.index
]



# 24. Create simplified gene column

# ChIPseeker reports the associated gene in geneId.
# The simplified 'gene' column is used by downstream scripts.

consensus.peaks.anno.df$gene <- consensus.peaks.anno.df$geneId



# 25. Move key columns to front


key.columns <- c(
  "peak_id",
  "gene",
  "annotation",
  "distanceToTSS",
  "avgAccess"
)


existing.key.columns <- key.columns[
  key.columns %in% colnames(consensus.peaks.anno.df)
]


other.columns <- setdiff(
  colnames(consensus.peaks.anno.df),
  existing.key.columns
)


consensus.peaks.anno.df <- consensus.peaks.anno.df[
  ,
  c(
    existing.key.columns,
    other.columns
  )
]



# 26. Save annotated consensus peak table


write.csv(
  consensus.peaks.anno.df,
  file = file.path(
    output.dir,
    "ATAC_annotated_consensus_peaks.csv"
  ),
  row.names = FALSE
)



# PART E: Identify peaks not returned by ChIPseeker




# 27. Identify consensus peaks without ChIPseeker annotation


unannotated.peak.id <- setdiff(
  consensus.peak.id,
  annotated.peak.id
)


cat(
  "\nConsensus peaks not returned by ChIPseeker:",
  length(unannotated.peak.id),
  "\n"
)



# 28. Save unannotated peak coordinates


unannotated.peaks <- consensus.peaks.table[
  consensus.peaks.table$peak_id %in% unannotated.peak.id,
]


write.csv(
  unannotated.peaks,
  file = file.path(
    output.dir,
    "ATAC_unannotated_consensus_peaks.csv"
  ),
  row.names = FALSE
)



# PART F: Summarize peak annotations

# 29. Count promoter-associated peaks


promoter.peaks <- grepl(
  "promoter",
  consensus.peaks.anno.df$annotation,
  ignore.case = TRUE
)


number.promoter.peaks <- sum(
  promoter.peaks,
  na.rm = TRUE
)


cat(
  "\nPromoter-associated consensus peaks:",
  number.promoter.peaks,
  "\n"
)

# Promoter-associated consensus peaks: 3973 

# 30. Count unique genes with accessible promoters


promoter.genes <- unique(
  consensus.peaks.anno.df$gene[
    promoter.peaks &
      !is.na(consensus.peaks.anno.df$gene)
  ]
)


cat(
  "Unique genes with at least one accessible promoter:",
  length(promoter.genes),
  "\n"
)

# Unique genes with at least one accessible promoter: 3568 

# 31. Calculate percent of annotated peaks in promoters


percent.promoter <- (
  100 *
    number.promoter.peaks /
    nrow(consensus.peaks.anno.df)
)


cat(
  "Percent of annotated consensus peaks classified as promoters:",
  round(percent.promoter, 2),
  "%\n"
)
# Percent of annotated consensus peaks classified as promoters: 35.81 %


# 32. Save annotation-category counts


annotation.counts <- as.data.frame(
  table(
    consensus.peaks.anno.df$annotation
  )
)


colnames(annotation.counts) <- c(
  "annotation",
  "number_of_peaks"
)


annotation.counts$percent_of_annotated_peaks <- (
  100 *
    annotation.counts$number_of_peaks /
    sum(annotation.counts$number_of_peaks)
)


annotation.counts <- annotation.counts[
  order(
    annotation.counts$number_of_peaks,
    decreasing = TRUE
  ),
]


write.csv(
  annotation.counts,
  file = file.path(
    output.dir,
    "ATAC_annotation_category_counts.csv"
  ),
  row.names = FALSE
)



# PART G: Final checks

# 33. Print final analysis summary


cat(
  "ATAC-seq samples:",
  nrow(samples),
  "\n"
)

cat(
  "Consensus definition: >= 2 of 4 samples\n"
)

cat(
  "MACS score threshold: 4.3\n"
)

cat(
  "Consensus peaks:",
  length(consensus.peaks),
  "\n"
)

cat(
  "Peaks returned by ChIPseeker:",
  nrow(consensus.peaks.anno.df),
  "\n"
)

cat(
  "Peaks not returned by ChIPseeker:",
  length(unannotated.peak.id),
  "\n"
)

cat(
  "Promoter-associated peaks:",
  number.promoter.peaks,
  "\n"
)

cat(
  "Unique genes with accessible promoters:",
  length(promoter.genes),
  "\n"
)

cat(
  "Percent promoter-associated:",
  round(percent.promoter, 2),
  "%\n"
)



# 34. Expected-result warnings


if (length(consensus.peaks) != 11231) {
  
  warning(
    "Expected 11,231 consensus peaks."
  )
  
}



if (number.promoter.peaks != 3973) {
  
  warning(
    "Expected approximately 3,973 promoter-associated consensus peaks."
  )
  
}



# 35. Save plain-text summary


capture.output(
  {
    
    cat(
      "ATAC-seq consensus peak construction\n"
    )
    

    
    cat(
      "ATAC-seq samples:",
      nrow(samples),
      "\n"
    )
    
    cat(
      "Consensus definition: peak present in at least 2 of 4 samples\n"
    )
    
    cat(
      "MACS -log10(p-value) threshold: 4.3\n"
    )
    
    cat(
      "Corresponding peak p-value threshold: 5 x 10^-5\n\n"
    )
    
    cat(
      "Consensus peaks:",
      length(consensus.peaks),
      "\n"
    )
    
    cat(
      "Peaks returned by ChIPseeker:",
      nrow(consensus.peaks.anno.df),
      "\n"
    )
    
    cat(
      "Peaks not returned by ChIPseeker:",
      length(unannotated.peak.id),
      "\n\n"
    )
    
    cat(
      "Promoter-associated peaks:",
      number.promoter.peaks,
      "\n"
    )
    
    cat(
      "Unique genes with accessible promoters:",
      length(promoter.genes),
      "\n"
    )
    
    cat(
      "Percent of annotated peaks classified as promoters:",
      round(percent.promoter, 2),
      "%\n"
    )
    
    cat(
      "\nMean CPM accessibility summary:\n"
    )
    
    print(
      summary(avg.accessibility)
    )
    
  },
  file = file.path(
    output.dir,
    "ATAC_consensus_peak_summary.txt"
  )
)



# 36. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "ATAC_consensus_peak_sessionInfo.txt"
  )
)
