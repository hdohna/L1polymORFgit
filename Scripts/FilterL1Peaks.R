##############################################
#
# General description:
#
#   The following script reads a L1 catalog, constructs genomic ranges for 
#   L1 and uses them to filter a bam file and filters reads that 
#   do overlap with reference L1s and writes out all the IDs 

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
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(ShortRead)
library(rtracklayer)
library(Rsamtools)
library(csaw)
library(GenomicRanges)

# Files and folders
L1CataloguePath <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
ChainFile       <- '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg38ToHg19.over.chain'
BamFilePath     <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19.bam'
FilteredBamFile <- gsub(".bam", ".filtered.bam", BamFilePath)

#######                       
# Get L1 ranges                    
#######                                     

cat("Getting reference L1 ranges \n")

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef          <- L1Catalogue$end_HG38 - L1Catalogue$start_HG38 > 5900
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1 & blnInRef),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
    ChainFilePath = ChainFile)
L1CatalogGR <- LiftOverList$GRCatalogue_hg19
L1CatalogGR <- resize(L1CatalogGR, 20000, fix = "center")

#####                                   
# Write IDs of reads in L1Ranges                        
#####

cat("Filtering bam file", BamFilePath, "\n")
paramFilter  <- ScanBamParam(mapqFilter = 1, which = L1CatalogGR)
filterBam(BamFilePath, FilteredBamFile, param = paramFilter)
