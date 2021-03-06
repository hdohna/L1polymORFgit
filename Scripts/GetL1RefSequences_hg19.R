##############################################
#
# General description:
#
#   The following script reads the repeat masker table 

# Input:
#
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
#source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)

# L1TableFileName   <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table_Hg19.csv"
# OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata'
L1TableFileName   <- "D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv"
ChrLPath          <- "D:/L1polymORF/Data/ChromLengthsHg19.Rdata"
OutResults        <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'

FullLength       <- 6000

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Loop through chromosomes and calculate length
# Chromosomes <- paste("chr", c(1:22, "X", "Y"), sep = "")
# ChromLengths <- sapply(Chromosomes, function(x) {
#   length(BSgenome.Hsapiens.UCSC.hg19[[x]])})
load(ChrLPath)
ChromLengths <- ChromLengthsHg19

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create names of subfamilies
repNChar    <- as.character(L1Table$repName)
SubFamiliesLookUp <- sapply(unique(repNChar), function(x){
  c(Name = x,
    SubFam = paste(strsplit(x, '[1-9]')[[1]][1:2], collapse = "1"))
})
NameMatch   <- match(repNChar, SubFamiliesLookUp["Name",])
SubFamilies <- SubFamiliesLookUp["SubFam", NameMatch]

#######################################
#                                     #
#    Create genomic ranges            #
#                                     #
#######################################

# Create GRanges objects with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)

# Get ranges of full-length L1HS
idxL1HsFullLength <- width(L1GRanges) >= FullLength & SubFamilies == "L1HS"
L1HSFullLength_GRanges <- L1GRanges[idxL1HsFullLength]

# Save results
cat("*******  Saving results ...   *******\n")
save(list = c("L1GRanges", "L1HSFullLength_GRanges", "ChromLengths"), 
     file = OutResults)

#######################################
#                                     #
#    Get sequences                    #
#                                     #
#######################################

# Get sequences of all LINE-1s
L1Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, L1GRanges)
