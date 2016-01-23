##############################################
#
# General description:
#
#   The following script reads a repeat table downloaded from the genome
#   browser repeatMasker track (http://genome.ucsc.edu/cgi-bin/hgTables)
#   and subsets to get all L1HS ranges

# Input:
#
#    D:/L1polymORF/Data/repeatsHg38: table with all repeats
#   

# Output:
#   
#    L1HS_repeat_table.csv: csv file with all L1HS repeats
#    L1Sequences_reptab.fas: fasta file with all L1 sequences

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Files and folders
RepeatFile          <- "D:/L1polymORF/Data/repeatsHg38"
TabOutfileName_L1HS <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"
SeqOutfileName_L1HS <- "D:/L1polymORF/Data/L1Sequences_reptab.fas"
TabOutfileName_L1   <- "D:/L1polymORF/Data/L1_repeat_table.csv"
SeqOutfileName_L1   <- "D:/L1polymORF/Data/L1AllSequences_reptab.fas"


#######################################
#                                     #
#     Read bed file and               #
#   save fasta file with sequences    #
#                                     #
#######################################

# Read repeat table
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]

######
# L1HS only
######

# Subset to get only L1HS rows and save
L1HSTable <- RepeatTable[RepeatTable$repName == "L1HS",]

# Subset to retain only strandard genome names
write.csv(L1HSTable, TabOutfileName_L1HS)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1HSIRanges <- IRanges(start = L1HSTable$genoStart,
                     end = L1HSTable$genoEnd)
L1HSGRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1HSIRanges,
                     strand = L1HSTable$strand)
L1HSGRanges <- L1HSGRanges[width(L1HSGRanges) >= 6000]

# Get all L1 sequences  
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1HSGRanges, as.character = T)
L1HSSeqNames <- paste(as.vector(seqnames(L1HSGRanges)), start(L1HSGRanges), end(L1HSGRanges),
      strand(L1HSGRanges), sep = "_")

# Write L1 sequences as fasta file
L1HSList <- lapply(L1HSSeq, s2c)
write.fasta(L1HSList, L1HSSeqNames, SeqOutfileName_L1HS)

######
# All L1s
######

# Subset to get only L1 rows and save
L1Table <- RepeatTable[RepeatTable$repFamily == "L1",]

# Subset to retain only strandard genome names
write.csv(L1Table, TabOutfileName_L1)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                       end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                       strand = L1Table$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000]

# Get all L1 sequences  
L1Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GRanges, as.character = T)
L1SeqNames <- paste(as.vector(seqnames(L1GRanges)), start(L1GRanges), end(L1GRanges),
                      strand(L1GRanges), sep = "_")

# Write L1 sequences as fasta file
L1List <- lapply(L1Seq, s2c)
write.fasta(L1List, L1SeqNames, SeqOutfileName_L1)

