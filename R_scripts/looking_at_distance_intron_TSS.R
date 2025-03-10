# Scratch: 1st exon length
library(rtracklayer)

library(dplyr)

# Read the GFF3 file
gff3 <- "~/Lab Notebook/Chapter2/Data_analysis/ATAC/Amil_ref_gen/Amil_v2.01/Amil.coding.gff3"
gff3.data <- read.table(gff3, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Assign column names to the GFF3 data
colnames(gff3.data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Extract relevant columns for exons
Amil.exons <- gff3.data %>% filter(type == "exon")

Amil.exons <- Amil.exons %>%
  mutate(gene.id = sub(".*Parent=([^;]+);.*", "\\1", attributes),
         exon.number = as.numeric(sub(".*exon:(\\d+);.*", "\\1", attributes)))

# Group by gene ID and select the exon with the lowest exon number for each gene
first.exons <- Amil.exons %>%
  group_by(gene.id) %>%
  filter(exon.number == min(exon.number)) %>%
  ungroup()

dim(first.exons) #28188

# Any duplicate genes?
dim(first.exons[duplicated(first.exons$gene.id), ]) # no

# Calculate the length of the first exons
first.exons$exon.length <- (first.exons$end - first.exons$start) + 1

# Is there a relationship between 1st exon length (proximity of 1st intron) and accessibility?
dim(all.genes.and.peaks)
all.genes.and.peaks$mean.norm.rc <- rowMeans(all.genes.and.peaks[, c("X1", "X2", "X3", "X4")], na.rm = TRUE)
all.genes.and.peaks$gene.id <- all.genes.and.peaks$transcriptId
all.genes.and.peaks <- all.genes.and.peaks %>% left_join(first.exons, by = join_by(gene.id))

all.genes.and.peaks <- all.genes.and.peaks %>%
  filter(!is.na(exon.length))

# Create the ggplot

lm_model <- lm(log10(mean.norm.rc) ~ log10(exon.length), data = all.genes.and.peaks)

# Get the slope and p-value from the linear model summary
lm_summary <- summary(lm_model)
slope <- lm_summary$coefficients[2, 1]
p_value <- lm_summary$coefficients[2, 4]

anova_result <- anova(lm_model)
anova_p_value <- anova_result["exon.length", "Pr(>F)"]


 ggplot(all.genes.and.peaks, aes(x = log10(exon.length), y = log10(mean.norm.rc), color = log10(baseMean))) +
  geom_point() +
   geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Exon Length", y = "Mean Normalized Read Count", color = "log10(Base Mean)") +
  theme_classic() +
   annotate("text", x = Inf, y = Inf, label = paste("Slope:", round(slope, 2), "\nANOVA P-value:", format.pval(anova_p_value, digits = 2)), 
            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
 

