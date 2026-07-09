#
#--------Lionshead soil fungi manuscript-----#
# Violet Harris / OSU / Kiser Lab / 2026
# 6.18.26

# This code rarefies the data post-DADA2.

# Reference materials: 
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

### CHECK FORMATTING ###
head(meta)
head(asv)

## Merge P and T reads within the same plots ##
# Remove the extra sample not in metadata/plots
counts <- counts[, colnames(counts) %in% meta$asv_name]
meta <- meta %>%
  filter(Plot != 0)

## Merge T and P within plots
merged_counts <- lapply(unique(meta$Plot), function(p) {
  samp_cols <- meta$asv_name[meta$Plot == p]
  rowSums(counts[, samp_cols, drop = FALSE])
})

## Convert to matrix/data frame
merged_asv <- as.data.frame(merged_counts)

## Use plot IDs as column names
names(merged_asv) <- unique(meta$Plot)

## Add ASV names back
merged_asv <- cbind(ASV = rownames(counts), merged_asv)

head(merged_asv)

# Matrix for rarefaction:
ASVs <- t(as.matrix(merged_asv[, -1]))

colnames(ASVs) <- merged_asv$ASV
rownames(ASVs) <- names(merged_asv)[-1]


# New plot-level metadata
meta_plot <- meta %>%
  group_by(Plot) %>%
  slice(1) %>%      # keep one row per plot
  ungroup() %>%
  select(-Sample, -PT, -asv_name)


# Save 
write.csv(meta_plot, "C:/Users/harrivio/LionsheadFungiData/inputs/mergedMETA.0626.csv")
saveRDS(ASVs, "C:/Users/harrivio/LionsheadFungiData/inputs/mergedASV.0626.rds")



### RAREFACTION ###
# Create function to rarefy asv table. Output is a table of each plot, asv, and abundance 
# after one rarefy iteration
rarefy <- function(x){
  rare_table <- rrarefy(x, raremax)
  rare_tibble <- (as_tibble(rare_table, rownames = NA) %>%
                    rownames_to_column(var = "Plot") %>%
                    pivot_longer(-Plot,
                                 names_to = "ASV",
                                 values_to = "Abundance"))
  return(rare_tibble)
}

# Test function and observe output
raremax <- min(rowSums(ASVs))
rarefy(ASVs)
raremax # The rarefaction threshold; 5165 in my case

# Iterate function to go through 1000 passes of rarefying the asv table
rarefy_iterations <- map_dfr(1:1000, ~rarefy(ASVs), .id = "iteration")

# Create table with mean values for each plot, asv, abundance combination
summary_iterations <- (rarefy_iterations %>%
                         group_by(Plot, ASV) %>%
                         summarize(mean_abund = mean(Abundance),
                                   .groups = "drop"))

##reformat the table so it becomes a ASV x plot matrix again
new_asv <- (summary_iterations %>%
              pivot_wider(names_from = ASV, values_from = mean_abund))
new_asv <- (new_asv %>%
              column_to_rownames(var = "Plot"))

# Create rarefaction curves to visualize the species accumulation.
rarecurvedata <- rarecurve(ASVs, step = 100)
map_dfr(rarecurvedata, bind_rows) %>%
  bind_cols(Group = rownames(ASVs),.) %>%
  pivot_longer(-Group) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name) %>%
  ggplot(aes(x=n_seqs, y=value, group = Group)) + geom_vline(xintercept = raremax, color = "gray") + geom_line()


## SAVE rarefied data frame
saveRDS(new_asv, "C:/Users/harrivio/LionsheadFungiData/inputs/rarefiedASVs.0626.rds")
