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
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(Rsamtools)

# Files and folders
BamFile            <- "/home/hzudohna/BoData/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
L1HSTableFileName  <- "/home/hzudohna/BoData/L1HS_repeat_table.csv"
OutputFileName     <- "/home/hzudohna/BoData/ReadsPerL1.RData"

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
                            
# Get reads for L1 ranges
param <- ScanBamParam(which = L1GRanges, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                        param = param)

# Save scanned reads as Rdata file
save.image(file = OutputFileName)

