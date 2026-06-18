#
#--------Lionshead soil fungi manuscript-----#
# Violet Harris / OSU / Kiser Lab / 2026
# 6.15.26

# This code rarefies the data post-DADA2.

# Code based on: 
# Abbey Neat: https://github.com/aneat/hja_fungi_manuscript/blob/main/code/01.remove_primers.R
# This tutorial: https://labolazar.github.io/amplicon_data_analysis/rarefaction.html#complete-code-1


# LOAD PACKAGES
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(dplyr)

# LOAD DATA
asv <- read.csv("C:/Users/harrivio/LionsheadFungiData/inputs/ASVtable.0626.csv")
tax <- read.csv("C:/Users/harrivio/LionsheadFungiData/inputs/TAXtable.0626.csv")
meta <- read.csv("C:/Users/harrivio/LionsheadFungiData/inputs/metadata_thesis_4.csv")

# Formatting
head(meta)
head(asv)

# Merge P and T reads within the same plots


