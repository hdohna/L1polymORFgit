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
BamFile   <- '/share/diskarray3/hzudohna/RNAseq/iPS0425-2_RNA-35180276.sorted.bam'
BamFolder <- '/share/diskarray3/hzudohna/RNAseq/'

# Specify file name parts
InputFileSuffix <- '.unique.sorted.bam'

# Get all bam files to be filtered
BamFilesToBeFiltered <- list.files(BamFolder, pattern = InputFileSuffix, 
                                  full.names = T)
BamFilesToBeFiltered <- BamFilesToBeFiltered[-grep(".bam.", 
                                                   BamFilesToBeFiltered)]

# Loop through bam files and filter them
for (BamFile in BamFilesToBeFiltered){
  cat("Filtering bam file", BamFile, "...\n")
  # WidthFilter <- FilterRules(minWidth <- function(DF){DF$qwidth > 100})
  # paramFilter  <- ScanBamParam(tagFilter = list(NM = 0:2, AS = 100:200),
  #                              mapqFilter = 1, what = c("mapq", "qwidth"))
  # filterBam(BamFile, FilteredBamFile, param = paramFilter, filter = WidthFilter)
  
  paramFilter  <- ScanBamParam(tagFilter = list(NM = 0:2, AS = 100:300),
                               mapqFilter = 1)
  FilteredBamFile <- gsub(".bam", ".filtered.bam", BamFile)
  filterBam(BamFile, FilteredBamFile, param = paramFilter)
  
}

