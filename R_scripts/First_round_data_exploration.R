## First round data exploration

# Download packages:
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)

# File Set-up
file.dir <- "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Gene_count_files/"
my.files <- grep(".txt", list.files(file.dir), value=TRUE)
my.metadata <- read.table("/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/tag_seq_meta_data.txt", header = TRUE)

# Create sample table
sampleNames <- c("07_88_Y", "03_72_R", "02_29_Y", "06_82_R", "01_51_Y", "08_95_R", "05_46_R", 
                 "05_43_R", "07_64_Y", "06_64_R", "03_48_Y", "07_39_R", "04_95_Y", "01_40_Y", 
                 "05_88_R", "01_52_Y", "01_45_Y", "02_34_Y", "08_50_R", "06_72_Y", "01_27_Y", 
                 "01_68_Y", "06_32_R", "08_30_R", "07_36_R", "03_29_R", "05_82_Y", "08_85_Y", 
                 "01_63_Y", "01_28_R", "07_53_Y", "04_65_Y", "02_65_R", "06_59_R", "04_48_R", 
                 "01_78_R", "05_50_Y", "04_70_Y", "05_62_Y", "08_33_Y", "01_46_Y", "01_54_Y", 
                 "08_45_R", "03_90_Y", "02_47_Y", "06_92_R", "02_79_R", "06_22_Y", "01_92_Y", 
                 "02_35_R", "01_42_R", "08_36_Y", "08_19_R", "03_91_R", "03_62_R", "01_58_Y",
                 "06_67_R", "07_31_R")

my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files, 
                             condition = my.metadata)



