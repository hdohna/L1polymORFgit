# The following script exports L1 ranges as bed file

# Load packages
library(ShortRead)
library(rtracklayer)
library(Rsamtools)
library(csaw)
library(GenomicRanges)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Read in table with known L1 
L1Catalogue <- read.csv('D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv', as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef          <- L1Catalogue$end_HG38 - L1Catalogue$start_HG38 > 5900
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1 & blnInRef),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
                                  ChainFilePath = 'D:/L1polymORF/Data/hg38ToHg19.over.chain')
L1CatalogGR <- LiftOverList$GRCatalogue_hg19
L1CatalogGRLarge <- resize(L1CatalogGR, 20000, fix = "center")

export.bed(L1CatalogGR, "D:/L1polymORF/Data/L1refRanges_hg19")