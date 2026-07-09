#
#--------Lionshead soil fungi manuscript-----#
# Violet Harris / OSU / Kiser Lab / 2026
# 6.15.26

# This code uses cutadapt to remove primers and dada2 to process ITS reads.

# Code based on: 
# Abbey Neat (adapted from Geoffry Zahn): https://github.com/aneat/hja_fungi_manuscript/blob/main/code/01.remove_primers.R
# Katie Hill (not available publicly -- sourced from Jed Cappellazzi in WSE at OSU CoF)
# DADA2 tutorial for ITS workflow: https://benjjneb.github.io/dada2/ITS_workflow.html


# LOAD PACKAGES
library(tidyverse)
library(purrr)
library(Biostrings)
library(ShortRead)
library(vegan)
library(reticulate)
library(dada2)
library(BiocParallel)
library(here)
library(ggplot2)

# the cutadapt tool is downloaded using python and miniconda, and both must be installed first if they are
# not already on the machine. Check online tutorials.



### PARSE FILE PATHS (name and sort forward/reverse reads) ####

path <- "C:/Users/harrivio/LionsheadFungiData/JustSequences" 
fnFs <- sort(list.files(path, pattern = "_R1_001\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001\\.fastq\\.gz$", full.names = TRUE))


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



### RUN CUTADAPT ###
cutadapt <- "C:/Users/harrivio/AppData/Local/Python/pythoncore-3.14-64/Scripts/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R
path.cut <- file.path(path, "cutadapt")

if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run cutadapt to remove the primers from your sequence
# for(i in seq_along(fnFs)) {
#   system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                              "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                              fnFs.filtN[i], fnRs.filtN[i])) # input files
# }

for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-m", "cutadapt",
                   R1.flags, R2.flags,
                   "-n", "2",
                   "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                   "-m", "1",
                   fnFs.filtN[i], fnRs.filtN[i]),
          stdout = TRUE, stderr = TRUE)
}
# Lots of output in the console is normal here while cutadapt is running.

# Sanity check: should return with 0 occurrences in any orientation. Check comes up with a few for my data, but not enough to worry about in the grand scheme.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



### Primer removal done. On to DADA2 for processing the raw ITS reads: ###


# Tell DADA2 which reads are forward and reverse
# This is to check the pattern, ours is diff than the documentation:
list.files(path.cut) 

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names) # Check that just the base name remains


# visualize a couple of read quality profiles to help select reasonable filtration parameters
plotQualityProfile(cutFs[1:2]) # forwards
plotQualityProfile(cutRs[1:2]) # reverse


### FILTER AND TRIM ###

# Filenames for filtered reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# This is the actual quality control step
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, #input and output names
                     maxN = 0, # requred for dada2 -- cannot have uncalled bases
                     maxEE = c(2, 4), # maximum expected errors allowed 
                     truncQ = 2, # standard cutoff for bad quality sequences
                     minLen = 50, # remove spurious low-length sequences
                     rm.phix = TRUE, # automatically remove PhiX spike-in reads from sequencing center ("control DNA" added during Illumina)
                     compress = TRUE, # compress output files with gzip
                     multithread = FALSE)  # on windows, set multithread = FALSE
head(out)

# save
saveRDS(out, "C:/Users/harrivio/LionsheadFungiData/dada2files/dada_out.0626.RDS")


### LEARN ERROR RATES ### 
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# save error files
saveRDS(errF, "C:/Users/harrivio/LionsheadFungiData/dada2files/errF.0626.RDS")
saveRDS(errR, "C:/Users/harrivio/LionsheadFungiData/dada2files/errR.0626.RDS")

# Sanity check: plot errors
errFplot <- plotErrors(errF, nominalQ = TRUE)
errRplot <- plotErrors(errR, nominalQ = TRUE)

# Save error plots
ggsave("C:/Users/harrivio/LionsheadFungiData/dada2files/errF_plots.png", 
       plot = errFplot, dpi = 500, height = 6, width = 6)
ggsave("C:/Users/harrivio/LionsheadFungiData/dada2files/errR_plots.png", 
       plot = errRplot, dpi = 500, height = 6, width = 6)


### DEREPLICATION  (De-replicate identical reads) ###
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#save dereplication files
saveRDS(derepFs, "C:/Users/harrivio/LionsheadFungiData/dada2files/derepF.0626.RDS")
saveRDS(derepRs, "C:/Users/harrivio/LionsheadFungiData/dada2files/derepR.0626.RDS")

### SAMPLE INFERENCE ###
set.seed(123)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# save sample inference files
saveRDS(dadaFs, "C:/Users/harrivio/LionsheadFungiData/dada2files/dadaFs.0626.RDS")
saveRDS(dadaRs, "C:/Users/harrivio/LionsheadFungiData/dada2files/dadaRs.0626.RDS")


### MERGE FWD AND REV READS ###
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# save mergers file
saveRDS(mergers, "C:/Users/harrivio/LionsheadFungiData/dada2files/mergers.0626.RDS")


### MAKE SEQUENCE TABLE (ASV table) ###
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

### REMOVE CHIMERAS ###
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "C:/Users/harrivio/LionsheadFungiData/dada2files/seqtab.nochim.0626.RDS")

# inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))


# SANITY CHECK: TRACK READS THROUGH PIPELINE to make sure everything worked
# there shouldn't be a step where most of the reads are removed
getN <- function(x) sum(getUniques(x))
out <- out[file.exists(filtFs),]  
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

saveRDS(track, "C:/Users/harrivio/LionsheadFungiData/dada2files/trackreads.0626.RDS")
write.csv(track, file = "C:/Users/harrivio/LionsheadFungiData/dada2files/trackreads.0626.csv")



### ASSIGN TAXONOMY ###
unite.ref <- "C:/Users/harrivio/DADA2_Thesis/UNITE_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# Check it out
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Save file
saveRDS(taxa, "C:/Users/harrivio/LionsheadFungiData/dada2files/taxa.0626.RDS")
write.csv(taxa, "C:/Users/harrivio/LionsheadFungiData/dada2files/taxa.0626.csv")


### MAKE FILES FOR LATER ###

# Defining objects and giving headers more manageable names w/ the "for" function
asv_seqs <- colnames(seqtab.nochim) #replace the column names
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


# Making files

#A. Make and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "C:/Users/harrivio/LionsheadFungiData/dada2files/asv_fasta.0626.fa")

#B. Make and write out a count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.csv(asv_tab, "C:/Users/harrivio/LionsheadFungiData/dada2files/ASVtable.0626.csv", quote=F)

#C. Make and write out a taxonomy table :
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.csv(asv_tax, "C:/Users/harrivio/LionsheadFungiData/dada2files/TAXtable.0626.csv", quote=F)


### END ###