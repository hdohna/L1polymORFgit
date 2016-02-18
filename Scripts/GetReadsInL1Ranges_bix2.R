
##############################################
#
# General description:
#
#   The following script reads in reads from a bam files that map to known
#   L1 ranges

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1HSTableFileName: path to file that contains L1HS ranges in a table


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
# source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(Rsamtools)
library(rtracklayer)

# Files and folders
BamFile            <- "/home/hzudohna/sorted_final_merged.bam"
L1HSTableFileName  <- "/home/hzudohna/L1polymORF/Data/L1HS_repeat_table.csv"
OutputFileName     <- "/home/hzudohna/L1polymORF/Data/ReadsPerL1.RData"
ChainFilePath      <- "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"

#######################################
#                                     #
#     Read L1 ranges and reads        #
#                                     #
#######################################

# Read in table with L1 ranges
L1Table <- read.csv(L1HSTableFileName)

# Create GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000] 

# Map to hg19 and retain only the coordinates that are uniquely mapped
L1GRanges19 <- liftOver(L1GRanges, chain = import.chain(ChainFilePath))
NrMapped    <- sapply(L1GRanges19, length)
idxMapped   <- which(NrMapped == 1)
L1GRanges19Mapped <- unlist(L1GRanges19[idxMapped])

# Get reads for L1 ranges
param <- ScanBamParam(which = L1GRanges19Mapped, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                        param = param)

# Save scanned reads as Rdata file
save.image(file = OutputFileName)

