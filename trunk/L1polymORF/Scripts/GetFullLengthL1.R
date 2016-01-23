##############################################
#
# General description:
#
#   The following script reads bed files of full-length L1s and the matching 
#   BSGenome to create a fasta file with all full-length L1 sequences

# Input:
#
#    hg38.fa_L1HS_6kb.bed: bed file with full-length L1 sequences
#   

# Output:
#   
#    L1Sequences.fas: fasta file with all L1 sequences

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Files and folders
DataFolder <- "D:/L1polymORF/Data/"
L1ReferenceBedFile <- "D:/L1polymORF/Data/hg38.fa_L1HS_6kb.bed"
RepeatFile<- "D:/L1polymORF/Data/repeatsHg38"
OutfileName <- "L1Sequences.fas"


#######################################
#                                     #
#     Read bed file and               #
#   save fasta file with sequences    #
#                                     #
#######################################

# Import BED files as GRanges      
cat("*******   Importing BED files as GRanges ...   *******\n")
L1Ref        <- import.bed(con = L1ReferenceBedFile) 

# Read repeat table
RepeatTable <- read.delim(RepeatFile)
AllNames <- unique(RepeatTable$repName)
grep("L1H", AllNames, value = T)
sum(RepeatTable$repName == "L1HS")

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Ref@elementMetadata$thick@width,
                     end = L1Ref@elementMetadata$thick@width +
                       L1Ref@elementMetadata$thick@start)
L1GRanges <- GRanges(seqnames = seqnames(L1Ref), ranges = L1IRanges,
                     strand = strand(L1Ref))

# Get all L1 sequences  
L1Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GRanges, as.character = T)
L1SeqNames <- paste(as.vector(seqnames(L1GRanges)), start(L1GRanges), end(L1GRanges),
      strand(L1GRanges), sep = "_")

# Write L1 sequences as fasta file
L1List <- lapply(L1Seq, s2c)
OutPath <- paste(DataFolder, OutfileName, sep = "")
write.fasta(L1List, L1SeqNames, OutPath)

