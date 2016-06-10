# The following script test for each L1 in the catalog whether it is present
# in a genome using PacBio data

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify flank width
Fwidth <- 500

# Boolean indicator whether reads from reference loci with high overlap should
# be filtered to a separate file:
blnFilterBam <- T

# Specify path to PacBio bam file
BamFilePath  <- "/share/diskarray3/hzudohna/sorted_final_merged.bam"
FilteredBamFilePath  <- "/share/diskarray3/hzudohna/PacBio/NA12878_filtered.bam"

# Read in table with known L1 
L1Catalog <- read.csv("/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                        as.is = T)
# L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv", 
#                       as.is = T)

# Get genomic ranges of L1 on hg19
blnMapped <- !is.na(L1Catalog$start_HG38) & L1Catalog$Allele == 1 
L1CatalogMapped <- L1Catalog[blnMapped, ]
# LiftOverList    <- LiftoverL1Catalog(L1CatalogMapped,
#                                      ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
LiftOverList    <- LiftoverL1Catalog(L1CatalogMapped)
L1CatalogMapped <- L1CatalogMapped[LiftOverList$idxUniqueMapped, ]
L1GRanges       <- LiftOverList$GRCatalogue_hg19
width(L1GRanges)


# Define flank ranges for getting reads intersecting with 
LeftFlankRanges  <- flank(L1GRanges, width = Fwidth, start = TRUE)
RightFlankRanges <- flank(L1GRanges, width = Fwidth, start = FALSE)

# Define parameters to scan barcodes from bam file
paramReadLeft  <- ScanBamParam(which = LeftFlankRanges, what = "qname")
paramReadRight <- ScanBamParam(which = RightFlankRanges, what = "qname")

# Scan reads from flanking regions
cat("Scan L1 flanking regions in", BamFilePath, "\n")
ReadIDListLeft  <- scanBam(BamFilePath, param = paramReadLeft)
ReadIDListRight <- scanBam(BamFilePath, param = paramReadRight)

# Loop through flanks and get the ratio of intersction to union
PropOverlap <- sapply(1:length(L1GRanges), function(i){
  LeftReadIDs <- ReadIDListLeft[[i]]$qname
  RightReadIDs <- ReadIDListRight[[i]]$qname
  length(intersect(LeftReadIDs, RightReadIDs)) / 
    length(union(LeftReadIDs, RightReadIDs))
})

# Plot a histogram of the proportion of overlapping reads 
pdf(file = "/home/hzudohna/L1polymORF/Figures/ProportionOverlappingReads_PacBio.pdf")
hist(PropOverlap, breaks = seq(0, 1, 0.05), 
     xlab = "Proportion of of overlapping reads",
     ylab = "Number of potential L1 loci",
     main = "")
dev.off()

# Add an indicator for L1 insertion
L1CatalogMapped$inNA12878 <- PropOverlap < 0.5
L1CatalogMapped$PropOverlap <- PropOverlap

# Filter out reads for igv inspection
if (blnFilterBam){
  cat("Filtering", BamFilePath, "\n")
  blnInReference <- (L1CatalogMapped$end_HG38 - 
                       L1CatalogMapped$start_HG38) > 6000
  blnHighOverlap <- PropOverlap > 0.5
  idx2Filter     <- which(blnInReference & blnHighOverlap)
  RangesToFilter <- c(LeftFlankRanges[idx2Filter], 
                      RightFlankRanges[idx2Filter])
  paramFilter    <- ScanBamParam(which = RangesToFilter)
  filterBam(BamFilePath, FilteredBamFilePath, param = paramFilter)
}

# Save results
save.image(file = "/home/hzudohna/L1polymORF/Data/L1PresenceTestResults_PacBio.RData")