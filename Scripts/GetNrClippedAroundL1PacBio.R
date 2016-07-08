# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify flank width
Fwidth <- 100

# Boolean indicator whether reads from reference loci with high overlap should
# be filtered to a separate file:
blnFilterBam <- T

# Specify path to PacBio bam file
BamFilePath      <- "/srv/gsfs0/projects/levinson/hzudohna/PacBio/sorted_final_merged.bam"
FilteredBamFilePath  <- "/srv/gsfs0/projects/levinson/hzudohna/PacBio/NA12878_filtered.bam"

# Read in table with known L1 
L1Catalog <- read.csv("/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                        as.is = T)

# Get genomic ranges of L1 on hg19
blnMapped <- !is.na(L1Catalog$start_HG38) & L1Catalog$Allele == 1 
L1CatalogMapped <- L1Catalog[blnMapped, ]
LiftOverList    <- LiftoverL1Catalog(L1CatalogMapped,
   ChainFilePath = "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg38ToHg19.over.chain")
L1CatalogMapped <- L1CatalogMapped[LiftOverList$idxUniqueMapped, ]
L1GRanges       <- LiftOverList$GRCatalogue_hg19
width(L1GRanges)

# Define flank ranges for getting reads intersecting with 
FlankRanges1  <- flank(L1GRanges, width = Fwidth, start = TRUE)
FlankRanges2  <- flank(L1GRanges, width = Fwidth, start = FALSE)
LeftFlankRanges  <- GRanges(seqnames = seqnames(LeftFlankRanges),
   ranges = IRanges(start = pmin(start(FlankRanges1), start(FlankRanges2)),
                    end = pmin(end(FlankRanges1), end(FlankRanges2))))
RightFlankRanges  <- GRanges(seqnames = seqnames(LeftFlankRanges),
   ranges = IRanges(start = pmax(start(FlankRanges1), start(FlankRanges2)),
                    end = pmax(end(FlankRanges1), end(FlankRanges2))))

# Create a vector of shift values
ShiftValues <- seq(0, 10*Fwidth, Fwidth)

# Loop over shift values and get the mean number of clipped bases per shift
ClipList <- lapply(ShiftValues, function(z){
  
  CurrentLeftRanges  <- shift(LeftFlankRanges, -z)
  CurrentRightRanges <- shift(RightFlankRanges, z)
  
  # Define parameters to scan barcodes from bam file
  paramReadLeft  <- ScanBamParam(which = CurrentLeftRanges, what = 'cigar')
  paramReadRight <- ScanBamParam(which = CurrentRightRanges, what = 'cigar')
  
  # Scan reads from flanking regions
  cat("Scan L1 flanking regions in", BamFilePath, "\n")
  ReadIDListLeft  <- scanBam(BamFilePath, param = paramReadLeft)
  ReadIDListRight <- scanBam(BamFilePath, param = paramReadRight)
  
  # Get the number of nucleotides that are clipped on the 5' end
  list(MeanClippedLeft = sapply(ReadIDListLeft, function(y){
    NrClipped <- sapply(y$cigar, function(x){
      CigNr          <- as.numeric(strsplit(x, "[A-Z]")[[1]])
      CigNr_5Prime   <- CigNr[1]
      CigNr_3Prime   <- CigNr[length(CigNr)]
      NrSplit        <- strsplit(x, "[1-9]")[[1]]
      CigType        <- grep("[A-Z]", NrSplit, value = T)
      CigType_5Prime <- CigType[1]
      CigType_3Prime <- CigType[max(c(1, length(CigType)))]
      c(CigNr_5Prime = CigNr_5Prime * (CigType_5Prime == "S"),
        CigNr_3Prime = CigNr_3Prime * (CigType_3Prime == "S"))
    }, USE.NAMES = F)
    rowMeans(NrClipped)
  }),
  MeanClippedRight = sapply(ReadIDListRight, function(y){
    NrClipped <- sapply(y$cigar, function(x){
      CigNr          <- as.numeric(strsplit(x, "[A-Z]")[[1]])
      CigNr_5Prime   <- CigNr[1]
      CigNr_3Prime   <- CigNr[length(CigNr)]
      NrSplit        <- strsplit(x, "[1-9]")[[1]]
      CigType        <- grep("[A-Z]", NrSplit, value = T)
      CigType_5Prime <- CigType[1]
      CigType_3Prime <- CigType[max(c(1, length(CigType)))]
      c(CigNr_5Prime = CigNr_5Prime * (CigType_5Prime == "S"),
        CigNr_3Prime = CigNr_3Prime * (CigType_3Prime == "S"))
    }, USE.NAMES = F)
    rowMeans(NrClipped)
  }))
})


# Save results
save.image(file = "/home/hzudohna/L1polymORFgit/Data/ClippedReadsL1_PacBio.RData")
