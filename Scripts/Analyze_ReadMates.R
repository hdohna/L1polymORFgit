##############################################
#
# General description:
#
#   The following script reads bam files of reads coming from chip-seq with
#   capture oligos containing L1HS sequences

# Input:
#
#    SRIP_eul1db: data on methods

# Output:
#   
#    : ...

##############################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(Rsamtools)

# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")
MinCoverage <- 5
IslRange    <- 1000
CalcNew_IslandPerCoverage <- F

#######################################
#                                     #
#    Turn BAM files into GRanges      #
#                                     #
#######################################

cat("*******   Turning BAM files into GRanges ...   *******\n")

