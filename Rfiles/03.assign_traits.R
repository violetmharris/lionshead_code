#
#--------Lionshead soil fungi manuscript-----#
# Violet Harris / OSU / Kiser Lab / 2026
# 7.9.26

# This annotates taxa with their functional guild traits from FUNGuild.
# This is post-rarefaction for my data (not necessary if read depths are relatively uniform).

# Reference materials: 
# FUNGuild: https://search.r-project.org/CRAN/refmans/MiscMetabar/html/funguild_assign.html
# Katie Hill (not available publicly -- sourced from Jed Cappellazzi in WSE at OSU CoF)

# LOAD PACKAGES
library(vegan)
library(dplyr)
library(tidyverse)
library(FUNGuildR)
library(stringr)
library(tibble)

# LOAD DATA
ASVs <- readRDS("C:/Users/harrivio/LionsheadFungiData/inputs/rarefiedASVs.0626.rds")
meta <- read.csv("C:/Users/harrivio/LionsheadFungiData/inputs/mergedMETA.0626.csv")
taxa <- read.csv("C:/Users/harrivio/LionsheadFungiData/dada2files/TAXtable.0626.csv")


# CHECK FORMATTING
colnames(meta)
colnames(ASVs)
rownames(ASVs)

# Fix funkiness in metadata
rownames(meta) <- meta$Plot # set plot numbers as rownames
meta <- meta[, -(1)] # remove unneeded column
meta <- meta[rownames(ASVs), ] # reorder rows to match ASV file

# Sanity checks
all(rownames(ASVs) %in% rownames(meta))
all(colnames(ASVs) %in% taxa$X)

# Save meta file to use in future
write.csv(meta, "C:/Users/harrivio/LionsheadFungiData/inputs/meta_formatted.0726.csv")



### PREP FUNGUILD-SPECIFIC FORMATTING REQUIREMENTS ###
# Vector of the rank columns b/c FUNGuild requires a single taxonomy string 
rank_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

# Clean NAs, make sure they are character rather than factor
taxa[rank_cols] <- lapply(taxa[rank_cols], function(x) {
  x <- as.character(x)
  x[x %in% c("<NA>", "NA", "")] <- NA
  x
})


# Replace NA with "" so paste won't insert "NA"
taxa[rank_cols] <- lapply(taxa[rank_cols], \(x) ifelse(is.na(x), "", x))

# Collapse multiple delimiters and trim trailing delimiters
collapse_taxonomy <- function(row_vals) {
  s <- paste(row_vals, collapse = ";")
  s <- gsub(";{2,}", ";", s)
  s <- sub("^;+","", s)
  s <- sub(";+$","", s)
  s
}

taxa$Taxonomy <- apply(taxa[, rank_cols, drop = FALSE], 1, collapse_taxonomy)

# check that it worked:
head(taxa)

# Create simplified input table:
funguild_input <- taxa[, c("X", "Taxonomy")]
colnames(funguild_input) <- c("ASV", "Taxonomy")

# Sanity check
head(funguild_input)



### FUNGuild Assignment ###

TraitTaxa <- funguild_assign(
  funguild_input,
  tax_col = "Taxonomy"
)

head(TraitTaxa) # sanity check

# Save
saveRDS(TraitTaxa, "C:/Users/harrivio/LionsheadFungiData/inputs/taxa_with_traits.0726.rds")
write.csv(TraitTaxa, "C:/Users/harrivio/LionsheadFungiData/inputs/taxa_with_traits.0726.csv")
