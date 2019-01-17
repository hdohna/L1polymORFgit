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
library(GenomicRanges)
library(rtracklayer)

ChrLPath    <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
OutResults  <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
OutBedPath  <- 'D:/L1polymORF/Data/L1HSRefRanges_hg19.bed'
OutBedPath2 <- 'D:/L1polymORF/Data/L1HSRefRanges_Plus200_hg19.bed'

FullLength       <- 6000

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Load vector with chromosome lengths
load(ChrLPath)

# Read in repeat masker table for hg1
RmskHg19 <- read.table("D:/L1polymORF/Data/repeats_Hg19", header = T, 
                       as.is = T, comment.char = "")

# Subset to get all L1
L1Table <- RmskHg19[RmskHg19$repFamily == "L1", ]

# Create a column with chromosome numbers
ChrChar <- as.character(L1Table$genoName)
L1Table$chromosome <- substr(ChrChar, 4, nchar(ChrChar))
table(nchar(L1Table$chromosome))
# Create genomic ranges for all L1s
L1GRanges_all <- makeGRangesFromDataFrame(L1Table,
                                      start.field = "genoStart", 
                                      end.field   = "genoEnd")

# Subset to get only L1HS or L1PA
blnL1HSPA <- substr(L1Table$repName, 1, 4) %in% c("L1HS", "L1PA")
blnChr <- nchar(L1Table$chromosome) <= 2
sum(blnL1HSPA)

# Export L1 ranges as bed files
export.bed(L1GRanges_all[blnChr], con = 'D:/L1polymORF/Data/L1allRefRanges_hg19.bed')
export.bed(L1GRanges_all[blnL1HSPA & blnChr], 
           con = 'D:/L1polymORF/Data/L1HSPARefRanges_hg19.bed')

# Write L1 table to csv file
write.csv(L1Table, 'D:/L1polymORF/Data/L1Table.csv')
