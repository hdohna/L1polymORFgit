# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and filters reads by various quality criteria

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify run parameters
blnCreateAlignList <- F

# Specify file paths
BamFile <- '/home/hzudohna/L1polymORF/Data/NA12878-L15P_S1_L001_001_Catalogue.dedup.unique.sorted.bam'
FilteredBamFile <- '/home/hzudohna/L1polymORF/Data/NA12878-L15P_aln2Catalogue_filtered.bam'
CatalogueFile <- "/home/hzudohna/L1polymORF/Data/L1CatalogUpdated_Fri_Apr_22_18-27-39_2016.csv"
AlignListFile <- "/home/hzudohna/L1polymORF/Data/CatlogueAlignList.RData"

############################
#                          #
#        Read Data         #
#                          #
############################

# Load L1 catalogue and files that indicate where in the alignment the L1 
# insertion is located
load(file = AlignListFile)

cat("Filtering bam file ...\n")
paramFilter  <- ScanBamParam(tagFilter = list(NM = 0:2),
                             scanBamFlag(isUnmappedQuery = F))
filterBam(BamFile, FilteredBamFile, param = paramFilter)
