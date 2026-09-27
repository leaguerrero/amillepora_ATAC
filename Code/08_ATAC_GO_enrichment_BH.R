# 08_ATAC_GO_enrichment_BH.R
#
# Purpose:
#   Test Gene Ontology enrichment among genes with
#   promoter-associated accessible chromatin regions in
#   Acropora millepora.
#
#   This script:
#     1. Loads the annotated ATAC-seq consensus peak set
#     2. Identifies genes with promoter-associated peaks
#     3. Defines the genome annotation background
#     4. Builds the gene-to-GO mapping
#     5. Performs classic Fisher GO enrichment separately for:
#          - Biological Process (BP)
#          - Molecular Function (MF)
#          - Cellular Component (CC)
#     6. Applies Benjamini-Hochberg correction
#     7. Retains terms with:
#          BH FDR < 0.05
#          at least 10 accessible-promoter genes in the term
#     8. Compares the BH-adjusted results with the original
#        nominal p < 0.01 criterion
#
# Inputs:
#   From Script 06:
#     - ATAC_annotated_consensus_peaks.csv
#
#   Genome annotation:
#     - Amil.all.maker.noseq.gff
#
#   Functional annotation:
#     - Amillepora_trinotate_annotation_report.csv
#
# Outputs:
#   - ATAC_Amil_GOmap.txt
#   - ATAC_GO_BP_all_terms.csv
#   - ATAC_GO_BP_BH_significant.csv
#   - ATAC_GO_MF_all_terms.csv
#   - ATAC_GO_MF_BH_significant.csv
#   - ATAC_GO_CC_all_terms.csv
#   - ATAC_GO_CC_BH_significant.csv
#   - ATAC_GO_BH_significant_all_ontologies.csv
#   - ATAC_GO_original_vs_BH_counts.csv
#   - ATAC_GO_terms_lost_after_BH.csv
#   - ATAC_GO_summary.txt
#
# Primary significance criterion:
#
#   BH FDR < 0.05
#   >= 10 accessible-promoter genes assigned to the GO term
#
# NOTE:
#   BH adjustment is performed separately within BP, MF,
#   and CC, matching the RNA-seq GO workflow.




# 1. Load packages
library(topGO)
library(GenomicFeatures)

# 2. Define input and output paths


project.dir <- "./Code_and_data/"


# Annotated consensus peaks from Script 06

annotated.peaks.file <- file.path(
  project.dir,
  "./Data/06_consensus_peaks",
  "ATAC_annotated_consensus_peaks.csv"
)


# A. millepora genome annotation

gff.file <- file.path(
  project.dir,
  "./external_resources/Amil_v2.01",
  "Amil.all.maker.noseq.gff"
)


# Trinotate functional annotation

trinotate.file <- file.path(
  project.dir,
  "./external_resources/Amil_v2.01",
  "Amillepora_trinotate_annotation_report.csv"
)


# Output directory

output.dir <- file.path(
  project.dir,
  "./Data/08_GO_enrichment"
)


dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# 3. Define GO significance criteria

fdr.threshold <- 0.05

minimum.accessible.genes <- 10

# PART A: Identify genes with accessible promoters

# 4. Read annotated consensus peaks

annotated.peaks <- read.csv(
  annotated.peaks.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Annotated consensus peak rows:",
  nrow(annotated.peaks),
  "\n"
)
# Annotated consensus peak rows: 11096 

# 5. Identify promoter-associated peaks

promoter.rows <- grepl(
  "promoter",
  annotated.peaks$annotation,
  ignore.case = TRUE
)


promoter.peaks <- annotated.peaks[
  promoter.rows,
]


cat(
  "Promoter-associated consensus peaks:",
  nrow(promoter.peaks),
  "\n"
)

# Promoter-associated consensus peaks: 3973 


# 6. Identify unique genes with accessible promoters


accessible.promoter.genes <- unique(
  promoter.peaks$gene[
    !is.na(promoter.peaks$gene) &
      promoter.peaks$gene != ""
  ]
)


# Remove transcript suffix if present

accessible.promoter.genes <- sub(
  "-RA$",
  "",
  accessible.promoter.genes
)


accessible.promoter.genes <- unique(
  accessible.promoter.genes
)


cat(
  "Unique genes with accessible promoters:",
  length(accessible.promoter.genes),
  "\n"
)

# Unique genes with accessible promoters: 3568 


# PART B: Define annotated genome background
# 7. Build transcript database from genome annotation


txdb <- makeTxDbFromGFF(
  gff.file,
  format = "gff3"
)



# 8. Extract annotated transcript IDs


a.millepora.transcripts <- transcripts(
  txdb
)


background.genes <- a.millepora.transcripts$tx_name


# Remove transcript suffix used by the genome annotation

background.genes <- sub(
  "-RA$",
  "",
  background.genes
)


# Retain unique non-missing identifiers

background.genes <- unique(
  background.genes[
    !is.na(background.genes) &
      background.genes != ""
  ]
)


cat(
  "\nAnnotated genome background genes:",
  length(background.genes),
  "\n"
)
# Annotated genome background genes: 38247 

# Original code produced 38,247 background identifiers

if (length(background.genes) != 38247) {
  
  warning(
    paste(
      "The previous analysis reported 38,247",
      "annotated background genes."
    )
  )
  
}

# 9. Confirm accessible-promoter genes belong to background

accessible.not.in.background <- setdiff(
  accessible.promoter.genes,
  background.genes
)


cat(
  "Accessible-promoter genes absent from genome background:",
  length(accessible.not.in.background),
  "\n"
)


# Restrict accessible-promoter set to the defined background

accessible.promoter.genes <- intersect(
  accessible.promoter.genes,
  background.genes
)


cat(
  "Accessible-promoter genes retained in background:",
  length(accessible.promoter.genes),
  "\n"
)

# Accessible-promoter genes retained in background: 3568 

# PART C: Build gene-to-GO mapping

# 10. Read Trinotate annotation

trinotate <- read.csv(
  trinotate.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "\nTrinotate annotation dimensions:",
  nrow(trinotate),
  "rows x",
  ncol(trinotate),
  "columns\n"
)

# 11. Define gene-ID column

colnames(trinotate)[1] <- "gene_id"


if (!"gene_ontology_blast" %in% colnames(trinotate)) {
  
  stop(
    paste(
      "Could not find gene_ontology_blast",
      "in the Trinotate annotation."
    )
  )
  
}

# 12. Standardize Trinotate gene IDs

trinotate$gene_id <- sub(
  "-RA$",
  "",
  trinotate$gene_id
)



# 13. Extract GO identifiers

#
# Trinotate stores GO terms together with annotation text.
# These substitutions retain the GO identifiers in the
# comma-separated format expected by topGO.


trinotate$GOlist <- gsub(
  "(GO:[0-9]*)\\^.*?\\`",
  "\\1,",
  trinotate$gene_ontology_blast
)


trinotate$GOlist <- gsub(
  "\\^.*?$",
  "",
  trinotate$GOlist
)


# 14. Build gene-to-GO table

GO.map <- trinotate[
  ,
  c(
    "gene_id",
    "GOlist"
  )
]


GO.map <- GO.map[
  !is.na(GO.map$gene_id) &
    GO.map$gene_id != "",
]


# 15. Save GO mapping

GO.map.file <- file.path(
  output.dir,
  "ATAC_Amil_GOmap.txt"
)


write.table(
  GO.map,
  file = GO.map.file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)


# 16. Read GO mapping into topGO

geneID2GO <- readMappings(
  file = GO.map.file,
  sep = "\t",
  IDsep = ","
)


cat(
  "Genes with GO mappings:",
  length(geneID2GO),
  "\n"
)


# 17. Examine GO annotation coverage

background.with.GO <- intersect(
  background.genes,
  names(geneID2GO)
)


accessible.with.GO <- intersect(
  accessible.promoter.genes,
  names(geneID2GO)
)


cat(
  "\nBackground genes represented in GO annotation:",
  length(background.with.GO),
  "\n"
)


cat(
  "Accessible-promoter genes represented in GO annotation:",
  length(accessible.with.GO),
  "\n"
)
# Accessible-promoter genes represented in GO annotation: 2766 


# PART D: Build topGO selection vector


# 18. Create named gene-selection vector

# 1 = gene has at least one promoter-associated consensus peak
# 0 = gene does not


accessible.gene.vector <- factor(
  as.numeric(
    background.genes %in% accessible.promoter.genes
  )
)


names(accessible.gene.vector) <- background.genes


cat(
  "\nAccessible-promoter gene coding:\n"
)


print(
  table(accessible.gene.vector)
)

# accessible.gene.vector
# 0     1 
# 34679  3568 

# PART E: Biological Process

# 19. Create Biological Process topGO object

GO.ATAC.BP <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = accessible.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)



# 20. Run classic Fisher test for Biological Process


BP.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


ATAC.Fisher.BP <- getSigGroups(
  GO.ATAC.BP,
  BP.test
)

# 21. Extract all Biological Process terms


ATAC.BP.all <- GenTable(
  GO.ATAC.BP,
  classic = ATAC.Fisher.BP,
  topNodes = length(
    ATAC.Fisher.BP@score
  ),
  numChar = 100
)

# 22. Convert Biological Process values to numeric


ATAC.BP.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    ATAC.BP.all$classic
  )
)


ATAC.BP.all$Significant <- as.numeric(
  ATAC.BP.all$Significant
)



# 23. Apply Benjamini-Hochberg correction


ATAC.BP.all$FDR <- p.adjust(
  ATAC.BP.all$classic_p,
  method = "BH"
)

# 24. Identify BH-significant Biological Process terms


ATAC.BP.significant <- ATAC.BP.all[
  !is.na(ATAC.BP.all$FDR) &
    ATAC.BP.all$FDR < fdr.threshold &
    ATAC.BP.all$Significant >= minimum.accessible.genes,
]


ATAC.BP.significant <- ATAC.BP.significant[
  order(ATAC.BP.significant$FDR),
]


cat(
  "\nBH-significant Biological Process terms:",
  nrow(ATAC.BP.significant),
  "\n"
)
# BH-significant Biological Process terms: 112 


# 25. Save Biological Process results


write.csv(
  ATAC.BP.all,
  file = file.path(
    output.dir,
    "ATAC_GO_BP_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  ATAC.BP.significant,
  file = file.path(
    output.dir,
    "ATAC_GO_BP_BH_significant.csv"
  ),
  row.names = FALSE
)



# PART F: Molecular Function
# 26. Create Molecular Function topGO object


GO.ATAC.MF <- new(
  "topGOdata",
  ontology = "MF",
  allGenes = accessible.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)

# 27. Run classic Fisher test for Molecular Function


MF.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


ATAC.Fisher.MF <- getSigGroups(
  GO.ATAC.MF,
  MF.test
)

# 28. Extract all Molecular Function terms


ATAC.MF.all <- GenTable(
  GO.ATAC.MF,
  classic = ATAC.Fisher.MF,
  topNodes = length(
    ATAC.Fisher.MF@score
  ),
  numChar = 100
)

# 29. Convert Molecular Function values to numeric


ATAC.MF.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    ATAC.MF.all$classic
  )
)


ATAC.MF.all$Significant <- as.numeric(
  ATAC.MF.all$Significant
)



# 30. Apply Benjamini-Hochberg correction


ATAC.MF.all$FDR <- p.adjust(
  ATAC.MF.all$classic_p,
  method = "BH"
)

# 31. Identify BH-significant Molecular Function terms


ATAC.MF.significant <- ATAC.MF.all[
  !is.na(ATAC.MF.all$FDR) &
    ATAC.MF.all$FDR < fdr.threshold &
    ATAC.MF.all$Significant >= minimum.accessible.genes,
]


ATAC.MF.significant <- ATAC.MF.significant[
  order(ATAC.MF.significant$FDR),
]


cat(
  "BH-significant Molecular Function terms:",
  nrow(ATAC.MF.significant),
  "\n"
)



# 32. Save Molecular Function results


write.csv(
  ATAC.MF.all,
  file = file.path(
    output.dir,
    "ATAC_GO_MF_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  ATAC.MF.significant,
  file = file.path(
    output.dir,
    "ATAC_GO_MF_BH_significant.csv"
  ),
  row.names = FALSE
)



# PART G: Cellular Component

# 33. Create Cellular Component topGO object


GO.ATAC.CC <- new(
  "topGOdata",
  ontology = "CC",
  allGenes = accessible.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)

# 34. Run classic Fisher test for Cellular Component


CC.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


ATAC.Fisher.CC <- getSigGroups(
  GO.ATAC.CC,
  CC.test
)


# 35. Extract all Cellular Component terms


ATAC.CC.all <- GenTable(
  GO.ATAC.CC,
  classic = ATAC.Fisher.CC,
  topNodes = length(
    ATAC.Fisher.CC@score
  ),
  numChar = 100
)



# 36. Convert Cellular Component values to numeric


ATAC.CC.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    ATAC.CC.all$classic
  )
)


ATAC.CC.all$Significant <- as.numeric(
  ATAC.CC.all$Significant
)



# 37. Apply Benjamini-Hochberg correction


ATAC.CC.all$FDR <- p.adjust(
  ATAC.CC.all$classic_p,
  method = "BH"
)



# 38. Identify BH-significant Cellular Component terms


ATAC.CC.significant <- ATAC.CC.all[
  !is.na(ATAC.CC.all$FDR) &
    ATAC.CC.all$FDR < fdr.threshold &
    ATAC.CC.all$Significant >= minimum.accessible.genes,
]


ATAC.CC.significant <- ATAC.CC.significant[
  order(ATAC.CC.significant$FDR),
]


cat(
  "BH-significant Cellular Component terms:",
  nrow(ATAC.CC.significant),
  "\n"
)
# BH-significant Cellular Component terms: 13 



# 39. Save Cellular Component results


write.csv(
  ATAC.CC.all,
  file = file.path(
    output.dir,
    "ATAC_GO_CC_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  ATAC.CC.significant,
  file = file.path(
    output.dir,
    "ATAC_GO_CC_BH_significant.csv"
  ),
  row.names = FALSE
)



# PART H: Combine BH-significant GO terms

# 40. Add ontology labels


ATAC.BP.significant$Type <- "BP"

ATAC.MF.significant$Type <- "MF"

ATAC.CC.significant$Type <- "CC"



# 41. Combine significant terms


ATAC.GO.significant <- rbind(
  ATAC.BP.significant,
  ATAC.MF.significant,
  ATAC.CC.significant
)


cat(
  "\nTotal BH-significant ATAC GO terms:",
  nrow(ATAC.GO.significant),
  "\n"
)

# 42. Put ontology column first


ATAC.GO.significant <- ATAC.GO.significant[
  ,
  c(
    "Type",
    setdiff(
      colnames(ATAC.GO.significant),
      "Type"
    )
  )
]



# 43. Save combined BH-significant table


write.csv(
  ATAC.GO.significant,
  file = file.path(
    output.dir,
    "ATAC_GO_BH_significant_all_ontologies.csv"
  ),
  row.names = FALSE
)




# Final summary
# 44. Print final analysis summary

cat(
  "Genome background genes:",
  length(background.genes),
  "\n"
)


cat(
  "Genes with accessible promoters:",
  length(accessible.promoter.genes),
  "\n"
)


cat(
  "Background genes with GO annotation:",
  length(background.with.GO),
  "\n"
)


cat(
  "Accessible-promoter genes with GO annotation:",
  length(accessible.with.GO),
  "\n\n"
)


cat(
  "BH FDR threshold:",
  fdr.threshold,
  "\n"
)


cat(
  "Minimum accessible-promoter genes per term:",
  minimum.accessible.genes,
  "\n\n"
)


cat(
  "BH-significant BP terms:",
  nrow(ATAC.BP.significant),
  "\n"
)


cat(
  "BH-significant MF terms:",
  nrow(ATAC.MF.significant),
  "\n"
)


cat(
  "BH-significant CC terms:",
  nrow(ATAC.CC.significant),
  "\n"
)


cat(
  "Total BH-significant GO terms:",
  nrow(ATAC.GO.significant),
  "\n"
)



# 45. Save plain-text summary


capture.output(
  {
    
    cat(
      "ATAC-seq promoter GO enrichment\n"
    )
    
    
    cat(
      "Background:\n"
    )
    
    cat(
      "All annotated A. millepora genes from genome annotation\n\n"
    )
    
    cat(
      "Genome background genes:",
      length(background.genes),
      "\n"
    )
    
    cat(
      "Genes with accessible promoters:",
      length(accessible.promoter.genes),
      "\n"
    )
    
    cat(
      "Background genes with GO annotation:",
      length(background.with.GO),
      "\n"
    )
    
    cat(
      "Accessible-promoter genes with GO annotation:",
      length(accessible.with.GO),
      "\n\n"
    )
    
    cat(
      "Primary significance criterion:\n"
    )
    
    cat(
      "Benjamini-Hochberg FDR <",
      fdr.threshold,
      "\n"
    )
    
    cat(
      "At least",
      minimum.accessible.genes,
      "accessible-promoter genes assigned to term\n\n"
    )
    
    cat(
      "BH-significant Biological Process terms:",
      nrow(ATAC.BP.significant),
      "\n"
    )
    
    cat(
      "BH-significant Molecular Function terms:",
      nrow(ATAC.MF.significant),
      "\n"
    )
    
    cat(
      "BH-significant Cellular Component terms:",
      nrow(ATAC.CC.significant),
      "\n"
    )
    
    cat(
      "Total BH-significant terms:",
      nrow(ATAC.GO.significant),
      "\n\n"
    )
    
    cat(
      "Comparison with original nominal p < 0.01 analysis:\n"
    )
    
    print(
      criterion.comparison
    )
    
    cat(
      "\nTerms lost after BH adjustment:",
      nrow(lost.after.BH),
      "\n"
    )
    
  },
  file = file.path(
    output.dir,
    "ATAC_GO_summary.txt"
  )
)



# 46. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "ATAC_GO_sessionInfo.txt"
  )
)
