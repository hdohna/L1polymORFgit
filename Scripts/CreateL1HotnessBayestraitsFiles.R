# The script below creates input files for a Bayestraits analysis. The files 
# are a nexus alignment of full-length L1 and a matching file of L1 hotness 
# values coded as a discrete character

################################
#                              #
#   Load packages and data     #
#                              #
################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)
library(ape)

# Load data about replication history (generated in script L1HS_repHistory.R)
load("D:/L1polymORF/Data/L1repHistoryResults.RData")

# Get table with L1 in reference genme from repeat masker
L1HSTable    <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table.csv", as.is = T)
L1HSTabNames <- paste(L1HSTable$genoName, L1HSTable$genoStart, L1HSTable$genoEnd,
                      L1Table$strand, sep = "_")

# Read prototypic L1 sequences (from http://www.girinst.org/repbase/)  
L1seqs <- read.dna("D:/L1polymORF/Data/Homo_sapiens_L1", format = "fasta",
                   as.character = T)

##########################################
#                                        #
#   Get genomic ranges and sequences     #
#                                        #
##########################################

# Make some corrections to create a proper GRanges object with L1 Seqences
L1HSIRanges <- IRanges(start = L1HSTable$genoStart,
                       end = L1HSTable$genoEnd)
L1HSGRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1HSIRanges,
                       strand = L1HSTable$strand)

# Subset ranges and table to get only full-length
subsetFullLength <- width(L1HSGRanges) >= 6000
L1HSGRanges      <- L1HSGRanges[subsetFullLength]
L1HSTable        <- L1HSTable[subsetFullLength, ]

# Get all L1 sequences  
cat("Get sequences for L1 element \n")
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1HSGRanges, as.character = T)
names(L1HSSeq) <- paste(as.vector(seqnames(L1HSGRanges)), start(L1HSGRanges), end(L1HSGRanges),
                      strand(L1HSGRanges), sep = "_")

# Get activity and total number of fragments
Act <- as.numeric(L1HSTable$Act_L1rp)

# Subset activity values and alignment to retain only those with measured
# activity
ActNotNA     <- Act[!is.na(Act)]
L1HSSeqNotNA <- L1HSSeq[!is.na(Act)]
L1HSSeqNotNA <- lapply(L1HSSeqNotNA, function(x){
  tolower(s2c(x))
})

##########################################
#                                        #
#   Add L1PA as root                     #
#                                        #
##########################################

# Remove \t from names of prototypic L1 sequences
names(L1seqs) <- gsub("\t", "_", names(L1seqs))
names(L1seqs) <- gsub(" ", "_", names(L1seqs))
SeqLengths    <- sapply(L1seqs, length)

# Find L1PA with maximum length
idxL1PA     <- grep("L1PA", names(L1seqs))
idxMaxL1PA  <- idxL1PA[which.max(SeqLengths[idxL1PA])]

# Add L1PA as root (first sequence)
L1HS_withRoot <- c(L1seqs[idxMaxL1PA], L1HSSeqNotNA)

# Write alignment as nexus file
write.fasta(L1HS_withRoot, names(L1HS_withRoot), 
            file.out = "D:/L1polymORF/Data/L1SequencesBayesTraitsUnaligned.fas")
run_MUSCLE(InputPath = "D:/L1polymORF/Data/L1SequencesBayesTraitsUnaligned.fas", 
           OutputPath = "D:/L1polymORF/Data/L1SequencesBayesTraitsAligned.fas")

# Read alignment and write it out as nexus file
L1HSAligned <- read.dna("D:/L1polymORF/Data/L1SequencesBayesTraitsAligned.fas",
                        format = "fasta", as.character = T)
write.nexus.data(L1HSAligned, "D:/L1polymORF/Data/L1SequencesBayesTraitsAligned.nex")

##########################################
#                                        #
#   Create Bayestraits input file        #
#                                        #
##########################################

# Add zero for the activitiy of L1PA
Act_withRoot <- c(0, ActNotNA)
names(Act_withRoot) <- names(L1HS_withRoot)

# Code values in three classes
Act_BayesTraits <- Act_withRoot
Act_BayesTraits[Act_withRoot == 0]                   <- 0L
Act_BayesTraits[Act_withRoot > 0 & Act_withRoot < 5] <- 1L
Act_BayesTraits[Act_withRoot >= 5]                   <- 2L

# Save Bayestraits file as text
write.table(Act_BayesTraits, file = "D:/L1polymORF/Data/ActivityBayesTraits.txt",
            col.names = F, quote = F)

##########################################
#                                        #
#   Create Bayestraits input file        #
#                                        #
##########################################

# Command for Bayestraits
# "C:\Program Files\BayesTraitsV2-Beta-Win64\BayesTraits.exe" D:\L1polymORF\Data\L1SequencesBayesTraitsAligned.nex.run1.t D:\L1polymORF\Data\ActivityBayesTraits.txt 

