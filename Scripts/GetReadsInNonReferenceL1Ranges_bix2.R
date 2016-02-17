
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

# Load results
load("/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata")

# Files and folders
BamFile            <- "/home/hzudohna/sorted_final_merged.bam"
OutputFileName     <- "/home/hzudohna/L1polymORF/Data/ReadsPerNonRefL1.RData"

#######################################
#                                                   #
#     Read non-reference L1 ranges and reads        #
#                                                   #
#######################################

# Extract indices for ranges from file names
idxRange <- sapply(FileNames, function(x) as.numeric(strsplit(x, "_")[[1]][2]))

# Create Ranges
L1NonRefRanges <- IslGRanges_reduced[idxRange]

# Match coordinates to hg19
L1NonRefRanges <- liftOver(L1NonRefRanges, 
   chain = import.chain("/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"))

# Get reads for L1 ranges
param <- ScanBamParam(which = L1NonRefRanges, what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                        param = param)

# Save scanned reads as Rdata file
save.image(file = OutputFileName)

