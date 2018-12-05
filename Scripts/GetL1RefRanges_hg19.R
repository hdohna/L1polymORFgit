##############################################
#
# General description:
#
#   The following script reads the repeat masker table for hg19
#   and creates genomic ranges of L1 and L1HS
#   with reference L1s 

# Input:
#
#     BamFile: path to file that contains mapped reads
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
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
# library(ShortRead)
# library(csaw)
# library(chipseq)
# library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(rtracklayer)

# L1TableFileName   <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table_Hg19.csv"
# OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata'
ChrLPath           <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
L1TableFileName <- "D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv"
OutResults      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
OutBedPath      <- 'D:/L1polymORF/Data/L1HSRefRanges_hg19.bed'
OutBedPath2     <- 'D:/L1polymORF/Data/L1HSRefRanges_Plus200_hg19.bed'

FullLength       <- 6000

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Load vector with chromosome lengths
load(ChrLPath)

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
save(list = c("L1GRanges", "L1HSFullLength_GRanges", "ChromLengthsHg19"), 
     file = OutResults)

# Export bed file with L1HS ranges
ChrChar <- as.character(L1Table$genoName)
ChrNr <- substr(ChrChar, 4, nchar(ChrChar))
L1GRanges_ChrNr <- GRanges(seqnames = ChrNr, ranges = L1IRanges,
                           strand = L1Table$strand)
export.bed(L1GRanges_ChrNr, con = OutBedPath)
L1GRangesPlus200 <- GRanges(seqnames = ChrNr, 
                            ranges = IRanges(start = start(L1GRanges) - 200,
                                      end = end(L1GRanges) + 200))
export.bed(L1GRangesPlus200, con = OutBedPath2)
