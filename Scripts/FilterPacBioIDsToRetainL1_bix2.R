##############################################
#
# General description:
#
#   The following script filters a PacBio bam file (NA12878) to retain only 
#   reads that map to L1

# Input:
#
#     BamFileToBeFiltered: path to bam file that has reads mapped to hg19
#     BamFilter: path to bam file that contains reads mapped to L1HS

# Output:
#   
#    : ...

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
BamFileToBeFiltered <- "/home/hzudohna/sorted_final_merged.bam"
BamFilter           <- "/home/hzudohna/NA12878PacBio_alnMappedSorted.bam"
OutFile             <- "/home/hzudohna/L1polymORF/Data/NA12878PacBio_L1hg19.bam"

# Files and folders (for testing purpose on laptop)
# BamFileToBeFiltered <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
# BamFilter           <- "D:/L1polymORF/Data/BamPerSuspectPeak/chr3_1161_aln.bam"
# OutFile             <- "D:/L1polymORF/Data/OutTest.bam"

#######                       
# Get L1 ranges                    
#######                                     

cat("Reading IDs of reads mapped to L1HS \n")

IDs <- scanBam(BamFilter, param=ScanBamParam(what="qname"))
IDs <- unlist(IDs)

cat("Filtering bam file by IDs of reads mapped to L1HS ...\n")
IDFilter <- FilterRules(getIDs <- function(DF){DF$qname %in% IDs})
filterBam(BamFileToBeFiltered, OutFile, filter=IDFilter)
