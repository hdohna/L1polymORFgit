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
Fwidth <- 1000

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
LeftFlankRanges  <- flank(L1GRanges, width = Fwidth, start = TRUE)
RightFlankRanges <- flank(L1GRanges, width = Fwidth, start = FALSE)

# Define parameters to scan barcodes from bam file
paramReadLeft  <- ScanBamParam(which = LeftFlankRanges, what = scanBamWhat())
paramReadRight <- ScanBamParam(which = RightFlankRanges, what = scanBamWhat())

# Scan reads from flanking regions
cat("Scan L1 flanking regions in", BamFilePath, "\n")
ReadIDListLeft  <- scanBam(BamFilePath, param = paramReadLeft)
ReadIDListRight <- scanBam(BamFilePath, param = paramReadRight)

# Get the number of nucleotides that are clipped on the 5' end
NrClipped <- sapply(ReadIDListLeft[[1]]$cigar, function(x){
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


# Loop through flanks and get the ratio of intersction to union
PropOverlap <- sapply(1:length(L1GRanges), function(i){
  LeftReadIDs <- ReadIDListLeft[[i]]$qname
  RightReadIDs <- ReadIDListRight[[i]]$qname
  length(intersect(LeftReadIDs, RightReadIDs)) / 
    length(union(LeftReadIDs, RightReadIDs))
})

# Calculate the proportion overlap based on read ranges
leftReads  <- lapply(LeftFlankRanges, function(x) extractReads(BamFilePath, x))
rightReads <- lapply(RightFlankRanges, function(x) extractReads(BamFilePath, x))
PropOverlapReadRanges <- sapply(1:length(leftReads), function(i){
  propLeft <- sum(overlapsAny(leftReads[[i]], RightFlankRanges[i])) / 
    length(leftReads[[i]])
  propRight <- sum(overlapsAny(rightReads[[i]], LeftFlankRanges[i])) / 
    length(rightReads[[i]])
  0.5 * (propLeft + propRight)
})
cor(PropOverlapReadRanges, PropOverlap)

# Plot a histogram of the proportion of overlapping reads 
pdf(file = "/home/hzudohna/L1polymORF/Figures/ProportionOverlappingReads_PacBio.pdf")
hist(PropOverlapReadRanges, breaks = seq(0, 1, 0.05), 
     xlab = "Proportion of of overlapping reads",
     ylab = "Number of potential L1 loci",
     main = "")
dev.off()

# Add an indicator for L1 insertion
L1CatalogMapped$inNA12878 <- PropOverlap < 0.5
L1CatalogMapped$PropOverlap <- PropOverlap

# Save results
save.image(file = "/home/hzudohna/L1polymORFgit/Data/L1PresenceTestResults_PacBio.RData")
