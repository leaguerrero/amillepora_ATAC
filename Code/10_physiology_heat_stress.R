# 10_physiology_heat_stress.R
#
# Purpose:
#   Analyze physiological responses to acute heat stress and
#   generate a 2 x 2 figure containing:
#
#     A. Fv/Fm
#     B. Red-channel intensity
#     C. Glutathione reductase activity
#     D. PCA of gene expression by acute heat-stress treatment
#
# Physiology panels:
#   - Day 0 is excluded because it preceded thermal-history
#     treatment.
#   - D07, D11, D14, and D21 are shown together.
#   - Each matched Control -> Heat Stress pair is connected.
#   - Connecting lines are colored by assay day.
#   - Point shape indicates thermal-history treatment.
#   - Black diamonds and error bars show model-adjusted means
#     and 95% confidence intervals.
#
# Final physiology model:
#
#   trait ~
#     Heat.Stress.Treatment +
#     Acc.Treatment +
#     Heat.Stress.Day +
#     (1 | Genotype)


# 1. Load packages

library(tidyverse)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)
library(DESeq2)
library(patchwork)

# 2. Define project and output directories

project.dir <- "./Code_and_data/"


output.dir <- file.path(
  project.dir,
  "./Data/10_heat_stress"
)


figure.output.dir <- file.path(
  project.dir,
  "Figures"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  figure.output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# PART A: PAM Fv/Fm
# 3. Read PAM/fragment metadata

pam.file <- file.path(
  project.dir,
  "./Data/00_inputs/physiology/",
  "PAM_FvFm_fragment_data.csv"
)


pam <- read.csv(
  pam.file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# 4. Prepare PAM data

pam <- pam %>%
  filter(
    !is.na(yield)
  ) %>%
  mutate(
    
    Genotype = factor(
      Genotype
    ),
    
    Acc.Treatment = factor(
      Acc.Treatment
    ),
    
    Heat.Stress.Day = factor(
      Heat.Stress.Day,
      levels = c(
        "D00",
        "D07",
        "D11",
        "D14",
        "D21"
      )
    ),
    
    Heat.Stress.Treatment = factor(
      Heat.Stress.Treatment,
      levels = c(
        "Control",
        "Heat Stress"
      )
    ),
    
    # PAM yield was stored x1000.
    # Convert to conventional Fv/Fm units.
    
    FvFm = yield / 1000,
    
    # Fragment identifier used to connect other assays
    # to the PAM metadata.
    
    FragKey = paste(
      as.integer(Genotype),
      as.integer(Tag.Number),
      toupper(
        substr(
          as.character(Color),
          1,
          1
        )
      ),
      sep = "_"
    ),
    
    # Matched Control / Heat Stress pair
    
    PairID = interaction(
      Genotype,
      Acc.Treatment,
      Heat.Stress.Day,
      drop = TRUE
    )
  )


cat(
  "PAM observations before removing Day 0:",
  nrow(pam),
  "\n"
)

# PART B: Red-channel intensity

# 5. Read RGB data

red.file <- file.path(
  project.dir,
  "./Data/00_inputs/physiology/",
  "red_channel_intensity_data.csv"
)


color.data <- read.csv(
  red.file,
  header = TRUE,
  stringsAsFactors = FALSE
)


cat(
  "\nRGB file columns:\n"
)


print(
  names(color.data)
)


# 6. Retain red-channel measurements

red.channel.data <- color.data %>%
  filter(
    !is.na(R)
  ) %>%
  mutate(
    red_intensity = R
  )


# 7. Parse red-channel fragment IDs

# Expected examples:
#
#   1-27Y
#   01_27_Y
#   01_27Y
#
# Result:
#
#   genotype_tag_color
#
# Example:
#
#   1_27_Y

red.fragment.label <- toupper(
  trimws(
    as.character(
      red.channel.data$Frag_id
    )
  )
)


red.fragment.parts <- str_match(
  red.fragment.label,
  "^(\\d+)[^0-9]*(\\d+)[^A-Z]*([RY])$"
)


red.channel.data$FragKey <- ifelse(
  is.na(red.fragment.parts[, 1]),
  NA_character_,
  paste(
    as.integer(
      red.fragment.parts[, 2]
    ),
    as.integer(
      red.fragment.parts[, 3]
    ),
    red.fragment.parts[, 4],
    sep = "_"
  )
)


cat(
  "\nRed-channel fragment IDs that could not be parsed:",
  sum(
    is.na(
      red.channel.data$FragKey
    )
  ),
  "\n"
)


if (any(is.na(red.channel.data$FragKey))) {
  
  print(
    red.channel.data[
      is.na(red.channel.data$FragKey),
    ]
  )
  
  stop(
    "Some red-channel fragment IDs could not be parsed."
  )
  
}


# 8. Create PAM metadata table for assay joins

pam.metadata <- pam %>%
  dplyr::select(
    FragKey,
    Genotype,
    Acc.Treatment,
    Heat.Stress.Day,
    Heat.Stress.Treatment,
    PairID
  )



# 9. Confirm PAM fragment keys are unique

pam.duplicate.keys <- pam.metadata %>%
  dplyr::count(
    FragKey
  ) %>%
  dplyr::filter(
    n > 1
  )

cat(
  "\nDuplicated PAM fragment keys:",
  nrow(pam.duplicate.keys),
  "\n"
)


if (nrow(pam.duplicate.keys) > 0) {
  
  print(
    pam.duplicate.keys
  )
  
  stop(
    "PAM fragment keys are not unique."
  )
  
}


# 10. Join red-channel measurements to experimental metadata

red.model.data <- red.channel.data %>%
  dplyr::select(
    FragKey,
    red_intensity
  ) %>%
  left_join(
    pam.metadata,
    by = "FragKey"
  )



# 11. Check red-channel metadata join


cat(
  "\nRed-channel observations:",
  nrow(red.model.data),
  "\n"
)

# Red-channel observations: 74 

cat(
  "Red-channel observations without matched PAM metadata:",
  sum(
    is.na(
      red.model.data$PairID
    )
  ),
  "\n"
)

# Red-channel observations without matched PAM metadata: 0 

if (any(is.na(red.model.data$PairID))) {
  
  print(
    red.model.data %>%
      filter(
        is.na(PairID)
      )
  )
  
  stop(
    "Some red-channel observations did not match PAM metadata."
  )
  
}


# 12. Restore factor structure

red.model.data <- red.model.data %>%
  mutate(
    
    Genotype = factor(
      Genotype
    ),
    
    Acc.Treatment = factor(
      Acc.Treatment
    ),
    
    Heat.Stress.Day = factor(
      Heat.Stress.Day,
      levels = c(
        "D00",
        "D07",
        "D11",
        "D14",
        "D21"
      )
    ),
    
    Heat.Stress.Treatment = factor(
      Heat.Stress.Treatment,
      levels = c(
        "Control",
        "Heat Stress"
      )
    )
  )


# PART C: Glutathione reductase

# 13. Read glutathione reductase data

gr.file <- file.path(
  project.dir,
  "./Data/00_inputs/physiology/",
  "GR_assay.csv"
)


gr.data <- read.csv(
  gr.file,
  header = TRUE,
  stringsAsFactors = FALSE
)


cat(
  "\nGR file columns:\n"
)


print(
  names(gr.data)
)


# 14. Parse GR fragment labels

# GR labels use a format such as:
#     1Y27
#
# which corresponds to:
#
#     genotype = 1
#     color    = Y
#     tag      = 27
#
# and therefore:
#
#     FragKey = 1_27_Y
#

gr.fragment.label <- toupper(
  trimws(
    as.character(
      gr.data$Label
    )
  )
)


gr.fragment.parts <- str_match(
  gr.fragment.label,
  "^(\\d+)([RY])(\\d+)$"
)


gr.data$FragKey <- ifelse(
  is.na(gr.fragment.parts[, 1]),
  NA_character_,
  paste(
    as.integer(
      gr.fragment.parts[, 2]
    ),
    as.integer(
      gr.fragment.parts[, 4]
    ),
    gr.fragment.parts[, 3],
    sep = "_"
  )
)


cat(
  "\nGR fragment IDs that could not be parsed:",
  sum(
    is.na(
      gr.data$FragKey
    )
  ),
  "\n"
)


if (any(is.na(gr.data$FragKey))) {
  
  print(
    gr.data[
      is.na(gr.data$FragKey),
    ]
  )
  
  stop(
    "Some GR fragment IDs could not be parsed."
  )
  
}


# 15. Join GR measurements to PAM metadata

gr.model.data <- gr.data %>%
  dplyr::select(
    FragKey,
    Glutathione_activity_normal
  ) %>%
  left_join(
    pam.metadata,
    by = "FragKey"
  )


# 16. Check GR metadata join

cat(
  "\nGR observations:",
  nrow(gr.model.data),
  "\n"
)

# GR observations: 74 

cat(
  "GR observations without matched PAM metadata:",
  sum(
    is.na(
      gr.model.data$PairID
    )
  ),
  "\n"
)

# GR observations without matched PAM metadata: 0 

if (any(is.na(gr.model.data$PairID))) {
  
  print(
    gr.model.data %>%
      filter(
        is.na(PairID)
      )
  )
  
  stop(
    "Some GR observations did not match PAM metadata."
  )
  
}


# 17. Prepare GR model data

gr.model.data <- gr.model.data %>%
  filter(
    !is.na(
      Glutathione_activity_normal
    )
  ) %>%
  mutate(
    
    Genotype = factor(
      Genotype
    ),
    
    Acc.Treatment = factor(
      Acc.Treatment
    ),
    
    Heat.Stress.Day = factor(
      Heat.Stress.Day,
      levels = c(
        "D00",
        "D07",
        "D11",
        "D14",
        "D21"
      )
    ),
    
    Heat.Stress.Treatment = factor(
      Heat.Stress.Treatment,
      levels = c(
        "Control",
        "Heat Stress"
      )
    )
  )


# 18. Log10-transform GR activity

if (
  any(
    gr.model.data$Glutathione_activity_normal <= 0,
    na.rm = TRUE
  )
) {
  
  gr.model.data$log10_GR_activity <- log10(
    gr.model.data$Glutathione_activity_normal +
      1e-6
  )
  
} else {
  
  gr.model.data$log10_GR_activity <- log10(
    gr.model.data$Glutathione_activity_normal
  )
  
}


# PART D: Remove Day 0
# 19. Exclude Day 0 from final physiology analysis

# Day 0 preceded the thermal-history treatment and therefore
# is not included in the final physiology models.


pam.post <- pam %>%
  filter(
    Heat.Stress.Day != "D00"
  ) %>%
  droplevels()


red.post <- red.model.data %>%
  filter(
    Heat.Stress.Day != "D00"
  ) %>%
  droplevels()


gr.post <- gr.model.data %>%
  filter(
    Heat.Stress.Day != "D00"
  ) %>%
  droplevels()


# 20. Check assay-day representation

cat(
  "\nPAM observations by assay day:\n"
)


print(
  table(
    pam.post$Heat.Stress.Day
  )
)
# D07 D11 D14 D21 
# 16  16  16  16 

cat(
  "\nRed-channel observations by assay day:\n"
)


print(
  table(
    red.post$Heat.Stress.Day
  )
)

# D07 D11 D14 D21 
# 16  16  16  16

cat(
  "\nGR observations by assay day:\n"
)


print(
  table(
    gr.post$Heat.Stress.Day
  )
)

# D07 D11 D14 D21 
# 16  16  16  16

# 21. Check Control / Heat Stress representation by day

cat(
  "\nPAM treatment by day:\n"
)


print(
  table(
    pam.post$Heat.Stress.Day,
    pam.post$Heat.Stress.Treatment
  )
)

# Control Heat Stress
# D07       8           8
# D11       8           8
# D14       8           8
# D21       8           8

cat(
  "\nRed treatment by day:\n"
)


print(
  table(
    red.post$Heat.Stress.Day,
    red.post$Heat.Stress.Treatment
  )
)

# Control Heat Stress
# D07       8           8
# D11       8           8
# D14       8           8
# D21       8           8

cat(
  "\nGR treatment by day:\n"
)


print(
  table(
    gr.post$Heat.Stress.Day,
    gr.post$Heat.Stress.Treatment
  )
)


# Control Heat Stress
# D07       8           8
# D11       8           8
# D14       8           8
# D21       8           8

# PART E: Final additive physiology models

# 22. PAM Fv/Fm model

pam.model <- lmer(
  FvFm ~
    Heat.Stress.Treatment +
    Acc.Treatment +
    Heat.Stress.Day +
    (1 | Genotype),
  data = pam.post,
  REML = FALSE
)


# 23. Red-channel intensity model

red.model <- lmer(
  red_intensity ~
    Heat.Stress.Treatment +
    Acc.Treatment +
    Heat.Stress.Day +
    (1 | Genotype),
  data = red.post,
  REML = FALSE
)


# 24. Glutathione reductase model

gr.model <- lmer(
  log10_GR_activity ~
    Heat.Stress.Treatment +
    Acc.Treatment +
    Heat.Stress.Day +
    (1 | Genotype),
  data = gr.post,
  REML = FALSE
)


# PART F: Type III model tests


# 25. PAM Type III ANOVA

pam.anova <- anova(
  pam.model,
  type = 3,
  ddf = "Satterthwaite"
)



cat(
  "PAM Fv/Fm\n"
)


print(
  summary(pam.model)
)


print(
  pam.anova
)


# 26. Red-channel Type III ANOVA

red.anova <- anova(
  red.model,
  type = 3,
  ddf = "Satterthwaite"
)



cat(
  "Red-channel intensity\n"
)


print(
  summary(red.model)
)


print(
  red.anova
)


# 27. GR Type III ANOVA

gr.anova <- anova(
  gr.model,
  type = 3,
  ddf = "Satterthwaite"
)


cat(
  "Glutathione reductase\n"
)


print(
  summary(gr.model)
)


print(
  gr.anova
)


# 28. Check model singularity

cat(
  "\nModel singularity checks:\n"
)


cat(
  "PAM:",
  isSingular(
    pam.model,
    tol = 1e-4
  ),
  "\n"
)


cat(
  "Red channel:",
  isSingular(
    red.model,
    tol = 1e-4
  ),
  "\n"
)


cat(
  "GR:",
  isSingular(
    gr.model,
    tol = 1e-4
  ),
  "\n"
)


# PART G: Model-adjusted acute heat-stress effects

# 29. PAM estimated marginal means

pam.stress.emm <- emmeans(
  pam.model,
  ~ Heat.Stress.Treatment,
  weights = "proportional"
)


pam.stress.means <- as.data.frame(
  confint(
    pam.stress.emm
  )
)


pam.stress.contrast <- as.data.frame(
  contrast(
    pam.stress.emm,
    method = list(
      "Heat Stress - Control" = c(
        -1,
        1
      )
    )
  )
)


cat(
  "\nPAM model-adjusted treatment means:\n"
)


print(
  pam.stress.means
)


cat(
  "\nPAM Heat Stress - Control:\n"
)


print(
  pam.stress.contrast
)


# 30. Red-channel estimated marginal means

red.stress.emm <- emmeans(
  red.model,
  ~ Heat.Stress.Treatment,
  weights = "proportional"
)


red.stress.means <- as.data.frame(
  confint(
    red.stress.emm
  )
)


red.stress.contrast <- as.data.frame(
  contrast(
    red.stress.emm,
    method = list(
      "Heat Stress - Control" = c(
        -1,
        1
      )
    )
  )
)


cat(
  "\nRed-channel model-adjusted treatment means:\n"
)


print(
  red.stress.means
)


cat(
  "\nRed-channel Heat Stress - Control:\n"
)


print(
  red.stress.contrast
)


# 31. GR estimated marginal means

gr.stress.emm <- emmeans(
  gr.model,
  ~ Heat.Stress.Treatment,
  weights = "proportional"
)


gr.stress.means <- as.data.frame(
  confint(
    gr.stress.emm
  )
)


gr.stress.contrast <- as.data.frame(
  contrast(
    gr.stress.emm,
    method = list(
      "Heat Stress - Control" = c(
        -1,
        1
      )
    )
  )
)


cat(
  "\nGR model-adjusted treatment means:\n"
)


print(
  gr.stress.means
)


cat(
  "\nGR Heat Stress - Control:\n"
)


print(
  gr.stress.contrast
)



# PART H: Save model results
# Save physiology model outputs


model.output.file <- file.path(
  output.dir,
  "physiology_model_outputs.txt"
)


capture.output(
  {
    
    # ========================================================
    # PAM Fv/Fm
    # ========================================================
    
    cat(
      "============================================================\n"
    )
    
    cat(
      "PAM Fv/Fm\n"
    )
    
    cat(
      "============================================================\n\n"
    )
    
    cat(
      "Model:\n"
    )
    
    cat(
      "FvFm ~ Heat.Stress.Treatment + Acc.Treatment + ",
      "Heat.Stress.Day + (1 | Genotype)\n\n"
    )
    
    
    cat(
      "MODEL SUMMARY\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      summary(pam.model)
    )
    
    
    cat(
      "\n\nTYPE III ANOVA\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      anova(
        pam.model,
        type = 3,
        ddf = "Satterthwaite"
      )
    )
    
    
    # ========================================================
    # Red-channel intensity
    # ========================================================
    
    cat(
      "\n\n============================================================\n"
    )
    
    cat(
      "Red-channel intensity\n"
    )
    
    cat(
      "============================================================\n\n"
    )
    
    cat(
      "Model:\n"
    )
    
    cat(
      "red_intensity ~ Heat.Stress.Treatment + Acc.Treatment + ",
      "Heat.Stress.Day + (1 | Genotype)\n\n"
    )
    
    
    cat(
      "MODEL SUMMARY\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      summary(red.model)
    )
    
    
    cat(
      "\n\nTYPE III ANOVA\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      anova(
        red.model,
        type = 3,
        ddf = "Satterthwaite"
      )
    )
    
    
    # ========================================================
    # Glutathione reductase
    # ========================================================
    
    cat(
      "\n\n============================================================\n"
    )
    
    cat(
      "Glutathione reductase activity\n"
    )
    
    cat(
      "============================================================\n\n"
    )
    
    cat(
      "Model:\n"
    )
    
    cat(
      "log10_GR_activity ~ Heat.Stress.Treatment + Acc.Treatment + ",
      "Heat.Stress.Day + (1 | Genotype)\n\n"
    )
    
    
    cat(
      "MODEL SUMMARY\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      summary(gr.model)
    )
    
    
    cat(
      "\n\nTYPE III ANOVA\n"
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(
      anova(
        gr.model,
        type = 3,
        ddf = "Satterthwaite"
      )
    )
    
  },
  
  file = model.output.file
)


cat(
  "\nPhysiology model outputs saved to:\n",
  model.output.file,
  "\n"
)

# PART I: Prepare physiology figure

# 34. Define assay-day colors
#
# Same colors are used in all three physiology panels.


day.colors <- c(
  "D07" = "#f94144",
  "D11" = "#f8961e",
  "D14" = "#90be6d",
  "D21" = "#577590"
)


day.labels <- c(
  "D07" = "Day 7",
  "D11" = "Day 11",
  "D14" = "Day 14",
  "D21" = "Day 21"
)


# 35. Shared physiology-plot theme

physiology.theme <- theme_classic(
  base_size = 19
) +
  theme(
    
    plot.title = element_text(
      size = 21,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 19,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16
    ),
    
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )


# PART J: Panel A — Fv/Fm

# 36. Plot Fv/Fm
# All assay days appear in the same panel.
#
# Line color = assay day
# Point shape = thermal-history treatment
#
# Black diamond = model-adjusted mean
# Black error bar = model-adjusted 95% CI

pam.plot <- ggplot(
  pam.post,
  aes(
    x = Heat.Stress.Treatment,
    y = FvFm,
    group = PairID
  )
) +
  
  geom_line(
    color = "grey85",
    alpha = 0.65,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      shape = Acc.Treatment
    ),
    color = "grey45",
    size = 3.4,
    alpha = 0.8
  ) +
  
  geom_errorbar(
    data = pam.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 2.0,
    color = "black"
  ) +
  
  geom_point(
    data = pam.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      y = emmean
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 6.0,
    stroke = 1.6,
    fill = "white",
    color = "black"
  ) +
  
  scale_color_manual(
    values = day.colors,
    breaks = names(day.colors),
    labels = day.labels,
    name = "Assay day"
  ) +
  
  labs(
    title = "Photosystem II efficiency",
    x = NULL,
    y = "Fv/Fm"
  ) +
  
  physiology.theme


pam.plot <- pam.plot +
  guides(shape = "none")

pam.plot


# PART K: Panel B — Red-channel intensity

# 37. Plot red-channel intensity

red.plot <- ggplot(
  red.post,
  aes(
    x = Heat.Stress.Treatment,
    y = red_intensity,
    group = PairID
  )
) +
  
  geom_line(
    color = "grey85",
    alpha = 0.65,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      shape = Acc.Treatment
    ),
    color = "grey45",
    size = 3.4,
    alpha = 0.8
  ) +
  
  geom_errorbar(
    data = red.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 2.0,
    color = "black"
  ) +
  
  geom_point(
    data = red.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      y = emmean
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 6.0,
    stroke = 1.6,
    fill = "white",
    color = "black"
  ) +
  
  labs(
    title = "Visual bleaching response",
    x = NULL,
    y = "Red-channel intensity"
  ) +
  
  physiology.theme

red.plot <- red.plot +
  guides(shape = "none")

red.plot


# PART L: Panel C — Glutathione reductase

# 38. Plot glutathione reductase
# Plot on the same log10 scale used for statistical inference.


gr.plot <- ggplot(
  gr.post,
  aes(
    x = Heat.Stress.Treatment,
    y = log10_GR_activity,
    group = PairID
  )
) +
  
  geom_line(
    color = "grey85",
    alpha = 0.65,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      shape = Acc.Treatment
    ),
    color = "grey45",
    size = 3.4,
    alpha = 0.8
  ) +
  
  geom_errorbar(
    data = gr.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 2.0,
    color = "black"
  ) +
  
  geom_point(
    data = gr.stress.means,
    aes(
      x = Heat.Stress.Treatment,
      y = emmean
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 6.0,
    stroke = 1.6,
    fill = "white",
    color = "black"
  ) +
  
  labs(
    title = "Glutathione reductase activity",
    x = NULL,
    y = expression(
      log[10] * " normalized GR activity"
    ),
    shape = "Thermal history"
  ) +
  
  physiology.theme

gr.plot <- gr.plot +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(
      t = -4,
      r = 0,
      b = 14,
      l = 0
    )
  )

gr.plot


# PART M: Panel D — whole-transcriptome PCA

# 39. Load final RNA-seq DESeq2 object

deseq.file <- file.path(
  project.dir,
  "./Data/03_DESeq2",
  "RNAseq_heat_stress_model.rds"
)


dds.heat <- readRDS(
  deseq.file
)


# 40. Perform variance-stabilizing transformation
#
# blind = TRUE because this PCA is descriptive QC /
# visualization rather than part of differential-expression
# model fitting.

vst.expression <- vst(
  dds.heat,
  blind = TRUE
)


# 41. Extract PCA data

pca.data <- plotPCA(
  vst.expression,
  intgroup = "heat_stress",
  returnData = TRUE
)


percent.variance <- round(
  100 *
    attr(
      pca.data,
      "percentVar"
    ),
  1
)


cat(
  "\nPCA variance explained:\n"
)


cat(
  "PC1:",
  percent.variance[1],
  "%\n"
)


cat(
  "PC2:",
  percent.variance[2],
  "%\n"
)


# 42. Set heat-stress treatment labels

pca.data$heat_stress <- factor(
  pca.data$heat_stress,
  levels = c(
    "24.5",
    "36"
  ),
  labels = c(
    "24.5°C control",
    "36°C heat stress"
  )
)


# 43. Plot gene-expression PCA

pca.plot <- ggplot(
  pca.data,
  aes(
    x = PC1,
    y = PC2,
    color = heat_stress
  )
) +
  
  geom_point(
    size = 4.3,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "24.5°C control" = "#bfd7ea",
      "36°C heat stress" = "#ff5a5f"
    ),
    name = "Heat-stress treatment"
  ) +
  
  labs(
    title = "Gene-expression PCA",
    x = paste0(
      "PC1 (",
      percent.variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent.variance[2],
      "%)"
    )
  ) +
  
  theme_classic(
    base_size = 19
  ) +
  
  theme(
    
    plot.title = element_text(
      size = 21,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 19,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16
    ),
    
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )


pca.plot


# PART N: Build 2 x 2 final figure

# 44. Combine panels


heat.response.figure <- (
  pam.plot |
    red.plot
) / (
  gr.plot |
    pca.plot
) +
  
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  
  theme(
    plot.tag = element_text(
      size = 22,
      face = "bold"
    ),
    legend.position = "bottom"
  )


heat.response.figure

# PART O: Save large publication-resolution figure

# 45. Save PDF

ggsave(
  filename = file.path(
    figure.output.dir,
    "heat_response_2x2_figure.pdf"
  ),
  plot = heat.response.figure,
  width = 16,
  height = 12,
  units = "in"
)


# 46. Save high-resolution TIFF

ggsave(
  filename = file.path(
    figure.output.dir,
    "heat_response_2x2_figure.tiff"
  ),
  plot = heat.response.figure,
  width = 16,
  height = 12,
  units = "in",
  dpi = 400,
  compression = "lzw"
)


# PART P: Save figure-input data

# 47. Save plotted physiology datasets


write.csv(
  pam.post,
  file = file.path(
    output.dir,
    "Figure_panel_A_FvFm_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  red.post,
  file = file.path(
    output.dir,
    "Figure_panel_B_red_channel_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  gr.post,
  file = file.path(
    output.dir,
    "Figure_panel_C_GR_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  pca.data,
  file = file.path(
    output.dir,
    "Figure_panel_D_gene_expression_PCA_data.csv"
  ),
  row.names = FALSE
)



# 48. Save session information

capture.output(
  sessionInfo(),
  file = file.path(
    output.dir,
    "physiology_sessionInfo.txt"
  )
)
