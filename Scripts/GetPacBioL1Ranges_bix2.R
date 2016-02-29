##############################################
#
# General description:
#
#   The following script reads in reads from PacBio bam file (NA12878) that
#   contains only reads mapped to L1 

# Input:
#
#     BamFile: path to bam file that has reads mapped to hg19 that were also 
#         mapped to the consensus L1HS sequence (NA12878PacBio_L1hg19.bam)

# Output:
#   
#    OutFile: GRanges of all reads aligned to L1 (PacBioL1Ranges.RData)

##############################################

######                                      
# Source packages and set parameters  
######

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
# source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(Rsamtools)
library(csaw)
library(GenomicRanges)

# Files and folders
BamFile  <- "/home/hzudohna/L1polymORF/Data/NA12878PacBio_L1hg19.bam"
OutFile  <- "/home/hzudohna/L1polymORF/Data/PacBioL1Ranges.RData"

# Load vector with chromosome lengths
load(file = "/home/hzudohna/L1polymORF/Data/ChromLengthsHg19.Rdata")

#######                       
# Get L1 ranges                    
#######                                     

cat("*******   Turning BAM file into GRanges ...   *******\n")

# Function to get read coverage per chromosome
ReadsPerChromPacBioL1 <- lapply(1:length(ChromLengthsHg19), function(i){
   Chrom       <- names(ChromLengthsHg19)[i]
   ChromLength <- ChromLengthsHg19[i]
   R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))
   cat("Extracting reads of chromosome", Chrom, "\n")
   Reads <- extractReads(bam.file = BamFile, region = R1)
})

# Save result
cat("*******   Saving results ...   *******\n")
save(list = c("ReadsPerChromPacBioL1"), file = OutFile)

