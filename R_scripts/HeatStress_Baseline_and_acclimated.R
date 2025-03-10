# Baseline heat stress

library(ggplot2)
library(ggpubr)
library(ggplot2)
library(plotly)
library(tidyverse)
library(stringr)
library(drc)
library(lmtest)
library(rstatix)
library(lme4)
#install.packages("Matrix")

# Plotting the difference in the PAM yields

# Load in the data
pam.yield <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/HeatStressExp/Acclimation Experiment Frag and Data Database - Coral Frags.csv")

pam.yield <-pam.yield[!is.na(pam.yield$yield), ]

my.metadata <- read.table("/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)
my.metadata$tube.label <- gsub("^0", "", my.metadata$ID)  # Remove leading 0
my.metadata$tube.label <- gsub("_", "-", my.metadata$tube.label)  # Replace underscore with hyphen
my.metadata$tube.label <- gsub("-(?=[A-Za-z])", "", my.metadata$tube.label, perl = TRUE)

# Baseline difference in thermal stress
## Subset day zero samples
D0.pam.yield <- subset(pam.yield, Heat.Stress.Day == "D00", 
                       select = c("Genotype", 
                                  "Heat.Stress.Day", 
                                  "Heat.Stress.Treatment", 
                                  "yield",
                                  "Tube.Label"))



# Plot the difference in yield between heat stress and control
compare_means(yield ~ Heat.Stress.Treatment, 
              data = D0.pam.yield, 
              method = "t.test",
              paired = TRUE)

baseline.heat.response.boxplot <- ggplot(D0.pam.yield, 
                                         aes(x = Heat.Stress.Treatment, 
                                             y = yield, 
                                             fill = Heat.Stress.Treatment)) +
                                  geom_boxplot() +
                                  labs(title = "Difference in PAM Yield between Heat Stress and Control Treatments",
                                       x = "Heat Stress Treatment",
                                       y = "Pam Yield") +
                                  scale_fill_manual(values = c("Control" = "#0a2463", 
                                                               "Heat Stress" = "#fb3640")) +
                                  theme_classic() +
                                  theme(legend.position = "none") +
                                  stat_compare_means(label = "p.signif", 
                                                     size = 6,
                                                     method = 't.test')

# ggsave("~/Lab Notebook/Chapter2/Data_analysis/HeatStressExp/Figures/baseline.hs.boxplot.jpeg", baseline.heat.response.boxplot)

# Model the proportion of PAM yield (HS/Control) for each day and have genotype 
# be a random effect in the model.
head(pam.yield)


pam.data <- na.omit(pam.yield[,c("Genotype",
                           "Tank",
                           "Acc.Treament",
                           "Heat.Stress.Day", 
                           "Heat.Stress.Treatment",
                           "yield",
                           "prop.yield",
                         "Tube.Label")])
pam.data$Genotype <- as.factor(pam.data$Genotype)
pam.data$Acc.class <- as.factor(pam.data$Acc.Treament)
pam.data$Heat.Stress.Day <- as.factor(pam.data$Heat.Stress.Day)
pam.data$acclimation.txt <- my.metadata[match(pam.data$Tube.Label, my.metadata$tube.label),4]
pam.data$acclimation.txt <- ifelse(pam.data$Tube.Label %in% c("5-63R", "5-90R", "6-38R", "6-60Y", "7-24R", "7-80Y", "8-20R", "8-57Y"), "control", 
                           ifelse(pam.data$Tube.Label %in% c("5-61R", "5-71R", "6-41R", "6-73Y", "7-43Y", "7-57R", "8-55Y", "8-94R"), "7DA_3DB", 
                                  pam.data$acclimation.txt))
# model the proportion yield~heat stress day with genotype as a random effect
# pam.yield.model <- lme4::lmer(prop.yield ~ Heat.Stress.Day + (1|Genotype), data = pam.data)
# summary(pam.yield.model)


pam.yield.model <- lm(prop.yield ~ Heat.Stress.Day * Acc.class, data = pam.data)
summary(pam.yield.model)
anova(pam.yield.model)

# Plot the relationship
acc.colors<- c("24.5" = "#696eff", "27" = "#e07a5f")

prop.pam.yield.scatterplot <- ggplot(pam.data, aes(x = Heat.Stress.Day, y = prop.yield, color = as.factor(Acc.Treament))) +
  geom_point() +
  geom_smooth(aes(x = as.numeric(Heat.Stress.Day), y = prop.yield),
              method = 'glm') +
  scale_color_manual(values = acc.colors, name = "Acclimation Treatment") +
  labs(title = "Scatter Plot of PAM yield proportion by heat stress day",
       x = "Heat Stress Day",
       y = "Proportion of PAM yield (HS/Control)") +
  theme_classic2()
#ggsave("~/Lab Notebook/Chapter2/Data_analysis/HeatStressExp/Figures/prop.yield.hs.day.scatterplot.jpeg",prop.pam.yield.scatterplot )



# RGB Color Analysis
## Code modified from Camille's archived code from Moorea
color.data <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/HeatStressExp/Acc-Rev Color Data - Sheet1.csv")

# find average value R, G, & B for each colony in each tank (across all replicates)
# average by group using dplyr
# average R
D0.color.data <- subset(color.data, Heat_stress_day== "D00")

# Plot the difference in R between heat stress and control
compare_means(R ~ Heat_stress_group, 
              data = D0.color.data, 
              method = "t.test",
              paired = TRUE)

baseline.heat.response.R.boxplot <- ggplot(D0.color.data, 
                                         aes(x = Heat_stress_group, 
                                             y = R, 
                                             fill = Heat_stress_group)) +
  geom_boxplot() +
  labs(title = "Difference in R channel intensity between Heat Stress and Control Treatments",
       x = "Heat Stress Treatment",
       y = "R Channel Intensity") +
  scale_fill_manual(values = c("control" = "#0a2463", 
                               "heat" = "#fb3640")) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif",
                     method = 't.test', 
                     paired = TRUE,
                     size = 6)

# Plot the difference in G between heat stress and control
compare_means(G ~ Heat_stress_group, 
              data = D0.color.data, 
              method = "t.test",
              paired = TRUE)

baseline.heat.response.G.boxplot <- ggplot(D0.color.data, 
                                           aes(x = Heat_stress_group, 
                                               y = G, 
                                               fill = Heat_stress_group)) +
  geom_boxplot() +
  labs(title = "Difference in G channel intensity between Heat Stress and Control Treatments",
       x = "Heat Stress Treatment",
       y = "G Channel Intensity") +
  scale_fill_manual(values = c("control" = "#0a2463", 
                               "heat" = "#fb3640")) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif",
                     method = 't.test', 
                     paired = TRUE,
                     size = 6)

# Plot the difference in B intensity between heat stress and control
compare_means(B ~ Heat_stress_group, 
              data = D0.color.data, 
              method = "t.test",
              paired = TRUE)

baseline.heat.response.B.boxplot <- ggplot(D0.color.data, 
                                           aes(x = Heat_stress_group, 
                                               y = B, 
                                               fill = Heat_stress_group)) +
  geom_boxplot() +
  labs(title = "Difference in B channel intensity between Heat Stress and Control Treatments",
       x = "Heat Stress Treatment",
       y = "B Channel Intensity") +
  scale_fill_manual(values = c("control" = "#0a2463", 
                               "heat" = "#fb3640")) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif",
                     method = 't.test', 
                     paired = TRUE,
                     size = 6)

# All three channel intensities in one plot
# Reshape the data frame for better plotting
reshaped.color.data <- color.data %>%
  gather(key = "Color", value = "Value", R, G, B)
reshaped.color.data$Acc_treatment <- as.factor(reshaped.color.data$Acc_treatment)
# Add p-values
reshaped.color.data.stats <- reshaped.color.data %>%
  group_by(Color) %>%
  t_test(Value ~ Heat_stress_group, paired = TRUE) %>% # Paired by genotype
  add_significance()
reshaped.color.data.stats
reshaped.color.data.stats <- reshaped.color.data.stats %>% add_xy_position(x = "Heat_stress_group")

# Create a box plot
color.intensity.bxp <- ggplot(reshaped.color.data, aes(x = Heat_stress_group, y = Value, fill = Color)) +
  geom_boxplot() +
  labs(title = "Box Plot of R, G, and B Intensity by Treatment",
       x = NULL,
       y = "Channel Intensity") +
  scale_fill_manual(values = c("R" = "red", "G" = "green", "B" = "blue")) +
  theme_classic()+
  theme(legend.position = "none") +
  facet_wrap(~Color) +
  stat_pvalue_manual(reshaped.color.data.stats)

#ggsave("~/Lab Notebook/Chapter2/Data_analysis/HeatStressExp/Figures/color.channel.intensity.bxp.jpeg", color.intensity.bxp)


# Add the proportion of the each color intensity of Heat:Control then plot the proportion

prop.R.scatterplot <- ggplot(reshaped.color.data, aes(x = Heat_stress_day, y = prop.R, color = Acc_treatment)) +
  geom_point() +
  geom_smooth(aes(x = as.numeric(Heat_stress_day), y = prop.R),
              method = 'glm') +
  scale_color_manual(values = acc.colors, name = "Acclimation Treatment") +
  labs(title = "Scatter Plot of R channel intensity proportion by heat stress day",
       x = "Heat Stress Day",
       y = "Proportion of R channel intensity (HS/Control)") +
  theme_classic2()

## Remove acclimation. Just heat stress vs control across all samples
pam.yield.hs.v.control <- subset(pam.yield, 
                       select = c("Genotype", 
                                  "Heat.Stress.Day", 
                                  "Heat.Stress.Treatment", 
                                  "yield",
                                  "Tube.Label"))

# Shapiro-Wilk test for normality
shapiro.test(pam.yield.hs.v.control$yield) # low p-value suggests data not normally dist.
qqnorm(pam.yield.hs.v.control$yield)
qqline(pam.yield.hs.v.control$yield)

compare_means(yield ~ Heat.Stress.Treatment, 
              data = pam.yield.hs.v.control, 
              method = "wilcox.test",
              paired = TRUE)

heat.response.boxplot <- ggplot(pam.yield.hs.v.control, 
                                         aes(x = Heat.Stress.Treatment, 
                                             y = yield, 
                                             fill = Heat.Stress.Treatment)) +
  geom_boxplot() +
  labs(title = "Fv/Fm yield by treatment",
       x = NULL,
       y = "Fv/Fm Yield") +
  scale_fill_manual(values = c("Control" = "#bfd7ea", 
                               "Heat Stress" = "#ff5a5f")) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", 
                     size = 6,
                     method = 'wilcox.test')  +
scale_x_discrete(labels = c("Control", "Heat")) 


heat.response.boxplot

## R channel intensity
red.chanel.data <- reshaped.color.data[reshaped.color.data$Color == "R",]

shapiro.test(red.chanel.data$Value) # high p-value suggests data normally dist.
qqnorm(red.chanel.data$Value)
qqline(red.chanel.data$Value)

# Paried T-test
control.r.chan.value <- red.chanel.data$Value[red.chanel.data$Heat_stress_group == "control"]
heat.r.chan.value <- red.chanel.data$Value[red.chanel.data$Heat_stress_group == "heat"]

t.test(control.r.chan.value, heat.r.chan.value, paired = TRUE)

R.chan.intensity.bxp <- ggplot(red.chanel.data, aes(x = Heat_stress_group, y = Value, fill = Heat_stress_group)) +
  geom_boxplot() +
  labs(title = "Visual bleaching score by treatment",
       x = "Heat Stress Treatment",
       y = "Visual Bleaching Score ") +
  scale_fill_manual(values = c("#bfd7ea", "#ff5a5f")) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", 
                     size = 6,
                     method = 't.test') +
  scale_x_discrete(labels = c("Control", "Heat"))

R.chan.intensity.bxp

((heat.response.boxplot /R.chan.intensity.bxp + plot_layout(guides = 'auto')) |  PCA_plot_hs) + plot_layout(guides = 'collect')

# Putting plots together:
patchwork <- ((heat.response.boxplot + R.chan.intensity.bxp)/PCA_plot_hs) 
patchwork + plot_annotation(
  title = expression(paste("Induced Bleaching Response of ", italic("Acropora millepora")))
)

# Glutathione reductase activity
gra.data <- read.csv("~/Lab Notebook/Chapter2/Data_analysis/Glutathione_reductase/GR_assay.csv", header = TRUE)
sample.meta <- read.table("~/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)

# Interested in Glutathione activity (normalized to total protein): 
## indicator for oxidative stress

# Preconditioning-based improvements in thermal tolerance 
# in Pocillopora acuta are accompanied by increases in 
# host glutathione reductase (GR) activity and
# gene expression (Majerova and Drury, 2022)

# We expect that GR activity in acclimated corals is higher in acclimated corals.
# and possibly decrease in the reverse-treatment corals. During heat stress, we 
# expect 

# We also expect the activity of GR activity to be influenced by heat-stress.
# From Majerova and Drury: LMM, activity~conditioning*time + (1| colony); p (time) = 0.001)

# GR Activity in acclimated fragments:
## Add acclimation treatment to the data frame:

gra.data$updated.label  <- str_replace(gra.data$Label, "^(\\d+)(\\D+)(\\d+)$", "0\\1_\\3_\\2")
gra.data$acclimation.txt <- sample.meta[match(gra.data$updated.label, sample.meta$ID),4]

# Note: Using the tag-seq meta data to pull acclimation treatment from, but we 
# didn't sequence a the 7DA_3DB samples, so entering that data manually.

gra.data$acclimation.txt <- ifelse(gra.data$updated.label %in% c("05_63_R", "05_90_R", "06_38_R", "06_60_Y", "07_24_R", "07_80_Y", "08_20_R", "08_57_Y"), "control", 
                                   ifelse(gra.data$updated.label %in% c("05_61_R", "05_71_R", "06_41_R", "06_73_Y", "07_43_Y", "07_57_R", "08_55_Y", "08_94_R"), "7DA_3DB", 
                                          gra.data$acclimation.txt))

custom.levels <- c("control", "7DA_0DB", "7DA_3DB", "7DA_7DB", "7DA_14DB")
gra.data$acclimation.txt <- factor(gra.data$acclimation.txt, levels = custom.levels)

gra.data <- gra.data %>%
  mutate(acclimation.group = ifelse(acclimation.txt == "control", "control", "acclimated"))

# Model the relationship
# glutathione reductase activity ~ heat.stress * acclimation + day
gra.model <- lm(Glutathione_activity_normal ~ Treatment * acclimation.txt + Day, data = gra.data)
summary(gra.model)
anova(gra.model)

shapiro.test(gra.data$Glutathione_activity_normal)
qqnorm(gra.data$Glutathione_activity_normal)
qqline(gra.data$Glutathione_activity_normal)

glut.red.boxplot <- ggplot(gra.data, aes(x = Treatment, y = Glutathione_activity_normal, fill = Treatment))+
  geom_boxplot() +
  scale_fill_manual(values = c("#bfd7ea", "#ff5a5f")) +
  labs(title = "Glutathione Reductase activity",
       x = "Heat Stress Treatment",
       y = "Glutathione Reductase Activity") +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", 
                     size = 6,
                     method = 'wilcox.test') +
  scale_x_discrete(labels = c("Control", "Heat")) 
glut.red.boxplot
# Putting plots together:
patchwork <- (heat.response.boxplot + R.chan.intensity.bxp + glut.red.boxplot) 
patchwork + plot_annotation(
  title = expression(paste("Induced Heat Response of ", italic("A. millepora")))
)

