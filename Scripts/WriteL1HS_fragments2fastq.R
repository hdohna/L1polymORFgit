##############################################
#
# General description:
#
#   The following script reads a repeat table downloaded from the genome
#   browser repeatMasker track (http://genome.ucsc.edu/cgi-bin/hgTables)
#   and writes all L1HS fragments as fastq file

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
library(ShortRead)
library(seqinr)

# Files and folders
RepeatFile   <- "D:/L1polymORF/Data/repeatsHg38"
FastqOutPath <- "D:/L1polymORF/Data/L1HS_Fragments.fastq"

# Minimum and maximum fragment size
MinFragSize <- 70
MaxFragSize <- 5000

#######################################
#                                     #
#     Read repeat table and               #
#   save fasta file with sequences    #
#                                     #
#######################################


# Read repeat table
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]

# Subset to get only L1HS rows with right fragment size
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS",]
genoRange   <- abs(RepeatTable$genoEnd - RepeatTable$genoStart)
RepeatTable <- RepeatTable[genoRange >= MinFragSize & genoRange <= MaxFragSize,]

# Make some corrections to create a proper GRanges object with L1 Seqences
L1HSIRanges <- IRanges(start = RepeatTable$genoStart,
                     end = RepeatTable$genoEnd)
L1HSGRanges <- GRanges(seqnames = RepeatTable$genoName, ranges = L1HSIRanges,
                     strand = RepeatTable$strand)

# Get sequences for L1 fragments 
cat("Get sequences for L1 element \n")
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1HSGRanges, as.character = T)
L1HSSeqNames <- paste(as.vector(seqnames(L1HSGRanges)), start(L1HSGRanges), end(L1HSGRanges),
      strand(L1HSGRanges), sep = "_")

# Write sequences out as fastq file
WriteFastq(L1HSSeq, L1HSSeqNames, FastqOutPath)
  


