# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and filters reads by various quality criteria

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
#source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify run parameters
blnCreateAlignList <- F

# Specify file paths
BamFile <- '/share/diskarray3/hzudohna/NA12878-L15P_S1_aln2Catalogue2016-05-07.dedup.unique.sorted.bam'
FilteredBamFile <- '/share/diskarray3/hzudohna/NA12878-L15P_aln2Catalogue2016-05-07_filtered.bam'
CatalogueFile <- "/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
AlignListFile <- "/home/hzudohna/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016_L1Locations.RData"

############################
#                          #
#        Read Data         #
#                          #
############################

# Load L1 catalogue and files that indicate where in the alignment the L1 
# insertion is located
#load(file = AlignListFile)

cat("Filtering bam file ...\n")
paramFilter  <- ScanBamParam(tagFilter = list(NM = 0:4),
                             mapqFilter = 1,
                             scanBamFlag(isUnmappedQuery = F))
filterBam(BamFile, FilteredBamFile, param = paramFilter)
