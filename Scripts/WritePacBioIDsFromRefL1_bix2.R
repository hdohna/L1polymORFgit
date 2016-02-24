##############################################
#
# General description:
#
#   The following script reads a bam file, finds reads that do overlap
#   with reference L1s and writes out all the IDs 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

######                                      
# Source packages and set parameters  
######

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(Rsamtools)
library(csaw)
library(GenomicRanges)

# Files and folders
L1TableFileName   <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table_Hg19.csv"
BamFile           <- "/home/hzudohna/NA12878PacBio_alnMapped.bam"
OutFilePath       <- "/home/hzudohna/L1polymORF/Data/PacBioIDsMapped2RefL1"

#######                       
# Get L1 ranges                    
#######                                     

cat("Getting reference L1 ranges \n")

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create GRanges objects with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)
L1GRanges[1]

#####                                   
# Write IDs of reads in L1Ranges                        
#####

cat("Writing reference L1 ranges \n")
WriteReadIDsInRanges(L1GRanges, InBamfilePath = BamFile, OutFilePath) 
  