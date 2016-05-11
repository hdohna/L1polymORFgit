# The following script reads a bam file of 10X reads aligned to L1, takes all
# the read IDs from the beginning and end of L1, gets the associated barcodes,
# constructs genomic ranges for barcodes and estimates L1 insertion sites


# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)
library(rtracklayer)

# Specify bam file paths for reads aligned to L1 and to the reference genome
L1BamPath  <- "/share/diskarray3/hzudohna/10XData/NA12878_10X_aln2L1.dedup.unique.sorted.bam"
RefBamPath <- "/share/diskarray3/hzudohna/10XData/NA12878_WGS_phased_possorted.bam"
FilterBamPath <- "/share/diskarray3/hzudohna/10XData/NA12878_10X_filterTemp.bam"

# Specify paths of reult files
BarCodeListFile <- "/share/diskarray3/hzudohna/10XData/Na12878_BarCodeListFile.RData"
BarCodeIDFile <- "/share/diskarray3/hzudohna/10XData/Na12878_BarcodesFullLength.RData"
BarCodeRangeFile <- "/share/diskarray3/hzudohna/10XData/Na12878_BarCodeRangeList.RData"

# Read in table with known L1 
L1Catalogue <- read.csv("/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv", 
                        as.is = T)

# Function to get all barcodes of reads that align to a range in L1
getBarcodes <- function(L1rangeStart, L1rangeEnd){
  
  # Define ranges to filter read IDs
  whichR <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
                 ranges = IRanges(start = L1rangeStart, end = L1rangeEnd))
  paramRead <- ScanBamParam(which = whichR, what = "qname")
  
  # Obtain read IDs of reads that intersect with specified range
  ReadIDs   <- scanBam(file = L1BamPath, param = paramRead)
  IDFilter    <- FilterRules(getIDs <- function(DF){DF$qname %in% ReadIDs})
  filterBam(RefBamPath, FilterBamPath, filter = IDFilter)
  paramReadFilter <- ScanBamParam(tag = "BX")# Need to specify ranges?
  BarcodeList <- scanBam(FilterBamPath, param = paramReadFilter)
  unique(unlist(BarcodeList))
}

######################
#                    #
#   Count barcodes   #
#                    #
######################

# # Define ranges to filter read IDs
# whichR <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
#                     ranges = IRanges(start = 1, end = 6000))
# paramRead <- ScanBamParam(which = whichR, what = "qname")
#   
# # Obtain read IDs of reads that intersect with specified range
# ReadIDs   <- scanBam(file = L1BamPath, param = paramRead)[[1]][[1]]
# IDFilter  <- FilterRules(getIDs <- function(DF){DF$qname %in% ReadIDs})
#   
# # Get barcodes of reads alingning to L1
# cat("***  Filtering reads aligning to L1 ...")
# filterBam(RefBamPath, FilterBamPath, filter = IDFilter)
# cat(" Done!   ***\n\n")
# cat("***  Getting barcodes of reads aligning to L1...")
# paramReadFilter <- ScanBamParam(tag = "BX")# Need to specify ranges?
# BarcodeList     <- scanBam(FilterBamPath, param = paramReadFilter)
# cat("***  Saving barcodes of reads alingning to bot L1 ends  \n")
# cat("in file", BarCodeListFile, "   ***\n\n")
# save(list = "BarcodeList", file = BarCodeListFile)
load(BarCodeListFile)

# # Get unique barcodes of reads aligning to 5' and 3' ends
# cat("***  Looking for barcodes of reads aligning to the L1 5' region ...")
# Barcodes5P <- getBarcodes(L1rangeStart = 1, L1rangeEnd = 1000)
# cat(" Done!   ***\n\n")
# cat("***  Looking for barcodes of reads aligning to the L1 3' region ...")
# Barcodes3P <- getBarcodes(L1rangeStart = 5000, L1rangeEnd = 6000)
# cat(" Done!   ***\n\n")
# BarcodesFullLength <- intersect(Barcodes5P, Barcodes3P)
# cat("***  Saving barcodes of reads alingning to bot L1 ends  \n")
# cat("in file", BarCodeIDFile, "   ***\n\n")
# save(list = "BarcodesFullLength", file = BarCodeIDFile)

# Loop through barcodes, scan reads per barcode, get minimum and maximum 
# positions among all reads from that barcode and put in into a GRanges
# object
cat("***  Determine genomic range per barcode   ***  \n")
Barcodes <- unlist(BarcodeList)
BarcodeTable <- table(Barcodes)
BarcodesFullLength <- names(BarcodeTable)[BarcodeTable > 1]
BarCodeRangeList <- lapply(BarcodesFullLength, function(x){
  
  paramFilter  <- ScanBamParam(tagFilter = list(BX = x),
                               what = c("rname", "pos"))# Need to specify ranges?
  ScannedReads <- scanBam(RefBamPath, FilterBamPath, paramFilter)
  idxStart <- which.min(ScannedReads[[1]]$pos)
  idxEnd   <- which.max(ScannedReads[[1]]$pos)
  Start    <- ScannedReads[[1]]$pos[idxStart]
  End      <- ScannedReads[[1]]$pos[idxEnd]
  
  GRanges(seqnames = ScannedReads[[1]]$rname[c(idxStart, idxEnd)],
          ranges = IRanges(start = Start, End = End))
})
cat("***  Saving barcode ranges in file", BarCodeRangeFile, "   ***\n")
save(list = "BarCodeRangeList", file = BarCodeRangeFile)

