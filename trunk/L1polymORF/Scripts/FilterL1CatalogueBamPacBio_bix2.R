# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and filters reads by various quality criteria

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify file paths
BamFile <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07sortedunique.sorted.bam'
FilteredBamFile <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07filtered.bam'

############################
#                          #
#        Read Data         #
#                          #
############################

cat("Filtering bam file ...\n")
paramFilter  <- ScanBamParam(tagFilter = list(AS = 5000:1000),
                             mapqFilter = 1,
                             scanBamFlag(isUnmappedQuery = F))
filterBam(BamFile, FilteredBamFile, param = paramFilter)
