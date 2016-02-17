
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
#    ReadsPerNonRefL1.RData: ...

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
cat("Load results from non-reference peak calling \n")
load("/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata")

# Files and folders
BamFile        <- "/home/hzudohna/sorted_final_merged.bam"
OutputFileName <- "/home/hzudohna/L1polymORF/Data/ReadsPerNonRefL1.RData"
ChainFilePath  <- "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"

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
L1NonRefRanges19 <- liftOver(L1NonRefRanges, 
                             chain = import.chain(ChainFilePath))
NrMapped <- sapply(L1NonRefRanges19, length)

# Retain only the coordinates that are uniquely mapped
L1NonRefRanges19Mapped <- unlist(L1NonRefRanges19[NrMapped == 1])
idxMapped <- idxRange[NrMapped == 1]

# Get reads for L1 ranges
cat("Read Pacbio reads per range \n")
param <- ScanBamParam(which = L1NonRefRanges19Mapped, 
                      what = scanBamWhat())
ScannedReads <- scanBam(file = BamFile, 
                        param = param)

# Save scanned reads as Rdata file
cat("Save results \n")
save.image(file = OutputFileName)

