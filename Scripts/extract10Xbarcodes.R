# The following script reads in reads from a 10X bam file and retains the barcodes
#  corresponding to reads aligning to the flanks of L1

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)
library(rtracklayer)

# Specify width of flanking regions to extract barcodes
Fwidth       <- 5000
FwidthFilter <- 10000

# Specify bam file path
InBamfilePath <- "/share/diskarray3/hzudohna/10XData/NA12878_WGS_phased_possorted.bam"
OutFilePrefix <- "/share/diskarray3/hzudohna/10XData/L1Flank5000_"
BamFolder     <- "/share/diskarray3/hzudohna/10XData/"
BamPrefix     <- "L1Flank5000Intersect_"
#OutFile <- "/share/diskarray3/hzudohna/10XData/test.bam"

# Read in table with known L1 
L1Catalogue <- read.csv("/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv", 
                        as.is = T)

# Get bam files with specified prefix, loop over files and extract accession numbers
BamFiles <- list.files(BamFolder, pattern = BamPrefix, full.names = T)
BamFiles <- BamFiles[grep(".bam", BamFiles)]
BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
AccProcessed <- sapply(BamFiles, function(x){
  Split <- strsplit(x, split = "_")[[1]][2]
  strsplit(Split, split = ".bam")[[1]][1]
})

# Create genomic ranges from catalogue
blnNotReference <- abs(L1Catalogue$end_HG38 - L1Catalogue$start_HG38) < 5000
#blnNotProcessed <- !L1Catalogue$Accession %in% AccProcessed
L1CatalogueMapped <- L1Catalogue[!is.na(L1Catalogue$start_HG38) &
                                   L1Catalogue$Allele == 1 &
                                   blnNotReference,]
GRCatalogue_hg38  <- GRanges(seqnames = L1CatalogueMapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogueMapped$start_HG38,
                                                     L1CatalogueMapped$end_HG38),
                                        end = pmax(L1CatalogueMapped$start_HG38,
                                                   L1CatalogueMapped$end_HG38)),
                       strand = L1CatalogueMapped$Strand)
GRCatalogue_hg19 <- liftOver(GRCatalogue_hg38, 
                      chain = import.chain(
                        "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"))
NrMapped_hg19    <- sapply(GRCatalogue_hg19, length)
idxUniqueMapped  <- which(NrMapped_hg19 == 1) 
GRCatalogue_hg19 <- unlist(GRCatalogue_hg38[idxUniqueMapped])

# Loop over L1s extract, barcodes per L1-flanking region, filter reads per barcode 
BarCodeList <- lapply(1:length(GRCatalogue_hg19), function(x){
  
  cat("Scanning barcodes for L1", x, "out of", length(GRCatalogue_hg19), "\n")
  
  # Get accession number
  AccNr <- L1CatalogueMapped$Accession[idxUniqueMapped[x]]
  
  # Define flank ranges for selecting barcodes 
  FlankRanges <- c(flank(GRCatalogue_hg19[x], width = Fwidth, start = TRUE),
                   flank(GRCatalogue_hg19[x], width = Fwidth, start = FALSE))
  
  # Define parameters to scan barcodes from bam file
  paramRead   <- ScanBamParam(which = FlankRanges, tag = "BX")
  
  # Scan barcodes from flanking regions
  BarCodesPerRangeList <- lapply(FlankRanges, function(y){
    paramRead   <- ScanBamParam(which = y, tag = "BX")
    scanBam(file = InBamfilePath, param = paramRead)
  })
  
  # Get unique barcodes that are not NA
  UniqueBarcodesLeft  <- unique(unlist(BarCodesPerRangeList[[1]]))
  UniqueBarcodesRight <- unique(unlist(BarCodesPerRangeList[[2]]))
  UniqueBarcodes      <- intersect(UniqueBarcodesLeft, UniqueBarcodesRight)
  UniqueBarcodes      <- UniqueBarcodes[!is.na(UniqueBarcodes)]
  
  # Set parameters to filter bam file
  if (length(UniqueBarcodes) > 0){
    paramFilter  <- ScanBamParam(tagFilter = list(BX = UniqueBarcodes))
    OutFile <- paste(OutFilePrefix, AccNr, ".bam", sep = "")
    filterBam(InBamfilePath, OutFile, param = paramFilter)
  }
})
