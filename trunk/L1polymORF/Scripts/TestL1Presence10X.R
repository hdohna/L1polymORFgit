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

# Specify path to PacBio bam file
BamFilePath  <- "/share/diskarray3/hzudohna/10XData/NA12878_WGS_phased_possorted.bam"

# Read in table with known L1 
L1Catalog <- read.csv("/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv",
                        as.is = T)

# Get genomic ranges of L1 on hg19
blnMapped <- !is.na(L1Catalog$start_HG38) & L1Catalog$Allele == 1 
L1CatalogMapped <- L1Catalog[blnMapped, ]
LiftOverList    <- LiftoverL1Catalog(L1CatalogMapped)
L1CatalogMapped <- L1CatalogMapped[LiftOverList$idxUniqueMapped, ]
L1GRanges       <- LiftOverList$GRCatalogue_hg19

# Define flank ranges for getting reads intersecting with 
LeftFlankRanges <- flank(L1GRanges, width = Fwidth, start = TRUE)
RightFlankRanges <- flank(L1GRanges, width = Fwidth, start = FALSE)

# Define parameters to scan barcodes from bam file
paramReadLeft  <- ScanBamParam(which = LeftFlankRanges, tag = "BX")
paramReadRight <- ScanBamParam(which = RightFlankRanges, tag = "BX")
# paramReadLeft  <- ScanBamParam(which = LeftFlankRanges, tag = "BC")
# paramReadRight <- ScanBamParam(which = RightFlankRanges, tag = "BC")

# Scan reads from flanking regions
ReadIDListLeft <- scanBam(BamFilePath, param = paramReadLeft)
ReadIDListRight <- scanBam(BamFilePath, param = paramReadRight)


# Loop through flanks and get the ratio of intersction to union
PropOverlap <- sapply(1:length(L1GRanges), function(i){
  LeftReadIDs <- ReadIDListLeft[[i]][[1]]$BX
  RightReadIDs <- ReadIDListRight[[i]][[1]]$BX
  length(intersect(LeftReadIDs, RightReadIDs)) / 
    length(union(LeftReadIDs, RightReadIDs))
})

# Plot a histogram of the proportion of overlapping reads 
pdf(file = "/home/hzudohna/L1polymORF/Figures/ProportionOverlappingReads_10X.pdf")
hist(PropOverlap, breaks = seq(0, 1, 0.05), 
     xlab = "Proportion of of overlapping barcodes",
     ylab = "Number of potential L1 loci",
     main = "")
dev.off()

# Add an indicator for L1 insertion
L1CatalogMapped$inNA12878 <- PropOverlap < 0.5
L1CatalogMapped$PropOverlap <- PropOverlap

# Save results
save.image(file = "/home/hzudohna/L1polymORF/Data/L1PresenceTestResults_10X.RData")
