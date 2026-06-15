#
#--------Lionshead soil fungi manuscript-----#
# Violet Harris / OSU / Kiser Lab / 2026
# 6.15.26

# This code uses cutadapt to remove primers. 

# Code based on: 
# Abbey Neat (adapted from Geoffry Zahn): https://github.com/aneat/hja_fungi_manuscript/blob/main/code/01.remove_primers.R
# Katie Hill (not available publicly -- sourced from Jed Cappellazzi in WSE at OSU CoF)
# DADA2 tutorial for ITS workflow: https://benjjneb.github.io/dada2/ITS_workflow.html


# LOAD PACKAGES
library(here)
library(tidyverse)
library(purrr)
library(Biostrings)
library(ShortRead)
library(vegan)
library(reticulate)
library(dada2)


### PARSE FILE PATHS (name and sort forward/reverse reads) ####

path <- "C:/Users/harrivio/DADA2_Thesis/JustSequences" 
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

# Sanity check 
print(fnFs)
print(fnRs)



### CHECK FOR AND REMOVE PRIMERS WITH CUTADAPT ###

# Supply primer sequences used during PCR
FWD <- "GTGAATCATCRAATYTTTG" # FWD primer sequence
REV <- "TCCTCCGCTTATTGATATGC" # REV primer sequence


# This function searches through both the forward and reverse reads and IDs where the primer input sequence is present
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  #The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Converts back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients

# Prefilter to remove reads with ambiguous (N) bases
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# Actually kick out the Ns. the maxN means that anything with over 0 Ns
# will be kicked out, the multithread has to do with processing power used
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)


# Finds primer matches, counts number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Checks the counts for the first N-filtered fastq file pair, you can see if 
# everything is properly aligned. See ITS workflow doc for comparison.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))



### Run cutadapt ###



