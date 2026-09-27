# 05_RNAseq_GO_enrichment.R
#
# Purpose:
#   Test Gene Ontology (GO) enrichment among genes that respond
#   to acute heat stress in Acropora millepora.
#
#   GO enrichment is performed separately for:
#     - Biological Process (BP)
#     - Molecular Function (MF)
#     - Cellular Component (CC)
#
#   Heat-responsive genes are compared against the background
#   of all genes that passed RNA-seq filtering and were tested
#   in the final DESeq2 model.
#
# Inputs:
#   From Script 03:
#     - RNAseq_heat_stress_model.rds
#     - RNAseq_heat_stress_DEGs.csv
#
#   Annotation:
#     - Amillepora_trinotate_annotation_report.csv
#
# Outputs:
#   - Amil_GOmap.txt
#   - GO_heat_BP_all_terms.csv
#   - GO_heat_BP_significant.csv
#   - GO_heat_MF_all_terms.csv
#   - GO_heat_MF_significant.csv
#   - GO_heat_CC_all_terms.csv
#   - GO_heat_CC_significant.csv
#   - GO_heat_significant_all_ontologies.csv
#   - GO_heat_FDR_threshold_comparison.csv
#   - GO_heat_summary.txt


# 1. Load packages
library(topGO)

# 2. Define input and output paths
project.dir <- "./Code_and_data/"


# DESeq2 model from Script 03

deseq.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_heat_stress_model.rds"
)


# Final heat-responsive gene set from Script 03

heat.genes.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_heat_stress_DEGs.csv"
)


# A. millepora Trinotate annotation

trinotate.file <- file.path(
  project.dir,
  "./external_resources/Amil_v2.01",
  "Amillepora_trinotate_annotation_report.csv"
)


# Output directory

output.dir <- file.path(
  project.dir,
  "./Data/05_GO_enrichment"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# 3. Define significance criteria

fdr.threshold <- 0.05

minimum.heat.genes <- 10

# 4. Load final DESeq2 model
dds.heat <- readRDS(
  deseq.file
)


background.genes <- rownames(
  dds.heat
)


cat(
  "Background genes from final DESeq2 analysis:",
  length(background.genes),
  "\n"
)

# Background genes from final DESeq2 analysis: 19649 

# Expected checkpoint

if (length(background.genes) != 19649) {
  
  warning(
    "Expected 19,649 genes in the RNA-seq background."
  )
  
}

# 5. Load final heat-responsive gene set

heat.results <- read.csv(
  heat.genes.file,
  stringsAsFactors = FALSE
)


heat.genes <- unique(
  heat.results$gene
)


cat(
  "Heat-responsive genes:",
  length(heat.genes),
  "\n"
)
# Heat-responsive genes: 2611 

# Expected checkpoint

if (length(heat.genes) != 2611) {
  
  warning(
    "Expected 2,611 final heat-responsive genes."
  )
  
}


# 6. Confirm that heat-responsive genes belong to background

heat.genes.not.in.background <- setdiff(
  heat.genes,
  background.genes
)


cat(
  "Heat-responsive genes absent from background:",
  length(heat.genes.not.in.background),
  "\n"
)


if (length(heat.genes.not.in.background) > 0) {
  
  warning(
    "Some heat-responsive genes are not present in the DESeq2 background."
  )
  
}

# PART A: Build gene-to-GO annotation map
# 7. Read A. millepora Trinotate annotation
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

# Trinotate annotation dimensions: 28188 rows x 16 columns

cat(
  "\nFirst several Trinotate column names:\n"
)

print(
  head(colnames(trinotate), 15)
)
# [1] "#gene_id"             "transcript_id"        "sprot_Top_BLASTX_hit"
# [4] "RNAMMER"              "prot_id"              "prot_coords"         
# [7] "sprot_Top_BLASTP_hit" "Pfam"                 "SignalP"             
# [10] "TmHMM"                "eggnog"               "Kegg"                
# [13] "gene_ontology_blast"  "gene_ontology_pfam"   "transcript"     


# 8. Define gene-ID column
# In the original Trinotate file, the first column contains
# the A. millepora gene identifier.

colnames(trinotate)[1] <- "gene_id"


# Confirm required GO column exists

if (!"gene_ontology_blast" %in% colnames(trinotate)) {
  
  stop(
    "Could not find the gene_ontology_blast column in the Trinotate annotation."
  )
  
}

# 9. Extract GO identifiers from Trinotate annotation
# The Trinotate GO field contains GO identifiers together with
# additional annotation text.
#
# These substitutions retain the GO IDs in the comma-separated
# format expected by topGO::readMappings().


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

# 10. Construct gene-to-GO mapping table


GO.map <- trinotate[
  ,
  c(
    "gene_id",
    "GOlist"
  )
]


# Remove rows without a gene ID

GO.map <- GO.map[
  !is.na(GO.map$gene_id) &
    GO.map$gene_id != "",
]


cat(
  "\nGenes in GO annotation table:",
  nrow(GO.map),
  "\n"
)

# Genes in GO annotation table: 28188 

# 11. Save GO map


GO.map.file <- file.path(
  output.dir,
  "Amil_GOmap.txt"
)


write.table(
  GO.map,
  file = GO.map.file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)



# 12. Read GO map into topGO

geneID2GO <- readMappings(
  file = GO.map.file,
  sep = "\t",
  IDsep = ","
)


cat(
  "Genes with GO mappings read by topGO:",
  length(geneID2GO),
  "\n"
)



# 13. Check annotation coverage of RNA-seq background


background.with.GO <- intersect(
  background.genes,
  names(geneID2GO)
)


heat.genes.with.GO <- intersect(
  heat.genes,
  names(geneID2GO)
)


cat(
  "\nRNA-seq background genes with GO annotation:",
  length(background.with.GO),
  "of",
  length(background.genes),
  "\n"
)
# RNA-seq background genes with GO annotation: 19434 of 19649 

cat(
  "Heat-responsive genes with GO annotation:",
  length(heat.genes.with.GO),
  "of",
  length(heat.genes),
  "\n"
)

# Heat-responsive genes with GO annotation: 2593 of 2611 

# PART B: Define heat-responsive genes for topGO
# 14. Create topGO gene-selection vector
# Every gene in the RNA-seq background receives:
#
#   1 = heat-responsive
#   0 = not heat-responsive
# ------------------------------------------------------------

heat.gene.vector <- factor(
  as.numeric(
    background.genes %in% heat.genes
  )
)


names(heat.gene.vector) <- background.genes


cat(
  "\nHeat-responsive gene coding for topGO:\n"
)

print(
  table(heat.gene.vector)
)
# heat.gene.vector
# 0     1 
# 17038  2611 


# PART C: Biological Process
# 15. Create Biological Process topGO object

GO.heat.BP <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = heat.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)



# 16. Run Fisher's exact test for Biological Process


BP.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


heat.Fisher.BP <- getSigGroups(
  GO.heat.BP,
  BP.test
)



# 17. Extract all Biological Process results


heat.BP.all <- GenTable(
  GO.heat.BP,
  classic = heat.Fisher.BP,
  topNodes = length(
    heat.Fisher.BP@score
  ),
  numChar = 100
)



# 18. Convert Biological Process p-values to numeric values
# topGO may display extremely small p-values as strings such as:
#
#     < 1e-30
#
# Remove the "<" before converting to numeric so these highly
# significant terms are not converted to NA.


heat.BP.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    heat.BP.all$classic
  )
)


# Benjamini-Hochberg multiple-testing correction

heat.BP.all$FDR <- p.adjust(
  heat.BP.all$classic_p,
  method = "BH"
)



# 19. Filter significant Biological Process terms


heat.BP.significant <- heat.BP.all[
  !is.na(heat.BP.all$FDR) &
    heat.BP.all$FDR < fdr.threshold &
    heat.BP.all$Significant >= minimum.heat.genes,
]


cat(
  "\nSignificant Biological Process terms:",
  nrow(heat.BP.significant),
  "\n"
)

# Significant Biological Process terms: 174 

# 20. Save Biological Process results


write.csv(
  heat.BP.all,
  file = file.path(
    output.dir,
    "GO_heat_BP_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  heat.BP.significant,
  file = file.path(
    output.dir,
    "GO_heat_BP_significant.csv"
  ),
  row.names = FALSE
)



# PART D: Molecular Function
# 21. Create Molecular Function topGO object


GO.heat.MF <- new(
  "topGOdata",
  ontology = "MF",
  allGenes = heat.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)



# 22. Run Fisher's exact test for Molecular Function


MF.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


heat.Fisher.MF <- getSigGroups(
  GO.heat.MF,
  MF.test
)



# 23. Extract all Molecular Function results


heat.MF.all <- GenTable(
  GO.heat.MF,
  classic = heat.Fisher.MF,
  topNodes = length(
    heat.Fisher.MF@score
  ),
  numChar = 100
)



# 24. Convert Molecular Function p-values to numeric values


heat.MF.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    heat.MF.all$classic
  )
)


# Benjamini-Hochberg multiple-testing correction

heat.MF.all$FDR <- p.adjust(
  heat.MF.all$classic_p,
  method = "BH"
)



# 25. Filter significant Molecular Function terms


heat.MF.significant <- heat.MF.all[
  !is.na(heat.MF.all$FDR) &
    heat.MF.all$FDR < fdr.threshold &
    heat.MF.all$Significant >= minimum.heat.genes,
]


cat(
  "Significant Molecular Function terms:",
  nrow(heat.MF.significant),
  "\n"
)
# Significant Molecular Function terms: 27 


# 26. Save Molecular Function results

write.csv(
  heat.MF.all,
  file = file.path(
    output.dir,
    "GO_heat_MF_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  heat.MF.significant,
  file = file.path(
    output.dir,
    "GO_heat_MF_significant.csv"
  ),
  row.names = FALSE
)



# PART E: Cellular Component


# 27. Create Cellular Component topGO object


GO.heat.CC <- new(
  "topGOdata",
  ontology = "CC",
  allGenes = heat.gene.vector,
  nodeSize = 10,
  annotationFun = annFUN.gene2GO,
  gene2GO = geneID2GO
)



# 28. Run Fisher's exact test for Cellular Component


CC.test <- new(
  "classicCount",
  testStatistic = GOFisherTest,
  name = "Fisher test",
  cutoff = 0.01
)


heat.Fisher.CC <- getSigGroups(
  GO.heat.CC,
  CC.test
)



# 29. Extract all Cellular Component results

heat.CC.all <- GenTable(
  GO.heat.CC,
  classic = heat.Fisher.CC,
  topNodes = length(
    heat.Fisher.CC@score
  ),
  numChar = 100
)



# 30. Convert Cellular Component p-values to numeric values


heat.CC.all$classic_p <- as.numeric(
  gsub(
    "^<\\s*",
    "",
    heat.CC.all$classic
  )
)


# Benjamini-Hochberg multiple-testing correction

heat.CC.all$FDR <- p.adjust(
  heat.CC.all$classic_p,
  method = "BH"
)



# 31. Filter significant Cellular Component terms


heat.CC.significant <- heat.CC.all[
  !is.na(heat.CC.all$FDR) &
    heat.CC.all$FDR < fdr.threshold &
    heat.CC.all$Significant >= minimum.heat.genes,
]


cat(
  "Significant Cellular Component terms:",
  nrow(heat.CC.significant),
  "\n"
)

Significant Cellular Component terms: 3 

# 32. Save Cellular Component results


write.csv(
  heat.CC.all,
  file = file.path(
    output.dir,
    "GO_heat_CC_all_terms.csv"
  ),
  row.names = FALSE
)


write.csv(
  heat.CC.significant,
  file = file.path(
    output.dir,
    "GO_heat_CC_significant.csv"
  ),
  row.names = FALSE
)



# PART F: Combine significant GO results




# 33. Add ontology labels


heat.BP.significant$Ontology <- "Biological Process"

heat.MF.significant$Ontology <- "Molecular Function"

heat.CC.significant$Ontology <- "Cellular Component"



# 34. Combine significant terms from all three ontologies

heat.GO.significant <- rbind(
  heat.BP.significant,
  heat.MF.significant,
  heat.CC.significant
)


cat(
  "\nTotal significant GO terms at FDR <",
  fdr.threshold,
  ":",
  nrow(heat.GO.significant),
  "\n"
)

# Total significant GO terms at FDR < 0.05 : 204 

# Save combined significant table

write.csv(
  heat.GO.significant,
  file = file.path(
    output.dir,
    "GO_heat_significant_all_ontologies.csv"
  ),
  row.names = FALSE
)


# PART G: Final summary
# 38. Print analysis summary


cat(
  "RNA-seq background genes:",
  length(background.genes),
  "\n"
)

cat(
  "Heat-responsive genes:",
  length(heat.genes),
  "\n"
)

cat(
  "Background genes with GO annotation:",
  length(background.with.GO),
  "\n"
)

cat(
  "Heat-responsive genes with GO annotation:",
  length(heat.genes.with.GO),
  "\n\n"
)

cat(
  "Using FDR <",
  fdr.threshold,
  "and at least",
  minimum.heat.genes,
  "heat-responsive genes per term:\n"
)

cat(
  "Biological Process:",
  nrow(heat.BP.significant),
  "\n"
)

cat(
  "Molecular Function:",
  nrow(heat.MF.significant),
  "\n"
)

cat(
  "Cellular Component:",
  nrow(heat.CC.significant),
  "\n"
)

cat(
  "Total:",
  nrow(heat.GO.significant),
  "\n"
)

# 39. Save plain-text analysis summary


capture.output(
  {
    
    cat(
      "Heat-responsive gene GO enrichment\n"
    )
    
    cat(
      "==================================\n\n"
    )
    
    cat(
      "RNA-seq background genes:",
      length(background.genes),
      "\n"
    )
    
    cat(
      "Heat-responsive genes:",
      length(heat.genes),
      "\n"
    )
    
    cat(
      "Background genes with GO annotation:",
      length(background.with.GO),
      "\n"
    )
    
    cat(
      "Heat-responsive genes with GO annotation:",
      length(heat.genes.with.GO),
      "\n\n"
    )
    
    cat(
      "Primary filtering threshold used in this script:\n"
    )
    
    cat(
      "BH FDR <",
      fdr.threshold,
      "\n"
    )
    
    cat(
      "At least",
      minimum.heat.genes,
      "heat-responsive genes assigned to term\n\n"
    )
    
    cat(
      "Significant terms:\n"
    )
    
    cat(
      "Biological Process:",
      nrow(heat.BP.significant),
      "\n"
    )
    
    cat(
      "Molecular Function:",
      nrow(heat.MF.significant),
      "\n"
    )
    
    cat(
      "Cellular Component:",
      nrow(heat.CC.significant),
      "\n"
    )
    
    cat(
      "Total:",
      nrow(heat.GO.significant),
      "\n\n"
    )
    
    cat(
      "Comparison of manuscript-relevant FDR thresholds:\n"
    )
    
    print(
      threshold.comparison
    )
    
  },
  file = file.path(
    output.dir,
    "GO_heat_summary.txt"
  )
)

# 40. Record R and package versions


capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "GO_heat_sessionInfo.txt"
  )
)
