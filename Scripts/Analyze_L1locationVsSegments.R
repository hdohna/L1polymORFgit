# The script below explore a catalog of full-length L1

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
#L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"

# Maximum fragment length
MaxFragLength <- 5900

######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

# Read in table with regulatory elements
SegmTable      <- read.table("D:/L1polymORF/Data/wgEncodeAwgSegmentationChromhmmGm12878", header = T)
#SegmTable      <- read.table("D:/L1polymORF/Data/wgEncodeAwgSegmentationCombinedChromhmmGm12878", header = T)
unique(SegmTable$name)
unique(SegmTable$itemRgb)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
SegmGR <- makeGRangesFromDataFrame(SegmTable, start.field = "ChromStart", 
                                      end.field = "ChromEnd")

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
L1RefGRFull <- L1RefGR[width(L1RefGR) > 6000]
L1FragmGR <- L1RefGR[width(L1RefGR) < MaxFragLength]
L1FragmNames <- paste(seqnames(L1FragmGR), start(L1FragmGR), end(L1FragmGR),
                      strand(L1FragmGR), sep = "_")

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create an overview table of L1 counts
L1CatLiftoverList <- LiftoverL1Catalog(L1CatalogL1Mapped, 
                                       ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR <- L1CatLiftoverList$GRCatalogue_hg19

##################################################
#                                                #
#    Calculate distances to genomic segments     #
#                                                #
##################################################

# Auxiliary function to get distances to closest gene
Dist2ClosestGR <- function(GR1, GR2){
  DistObj <- distanceToNearest(GR1, GR2, ignore.strand = T) 
  DistObj@elementMetadata@listData$distance
}


# Calculate distances from full-length L1 to nearest feature
L1Dist_cat <- Dist2ClosestGR(L1CatalogGR, SegmGR)
L1Dist_fragm <- Dist2ClosestGR(L1FragmGR, SegmGR)
sum(L1Dist_fragm == 0)


# Test for enrichment in different segment classes
EnrichTests <- Compare2RangesWith3rd(L1CatalogGR, L1FragmGR, SegmGR, 
                                     SegmTable$name)
barplot(EnrichTests$ORs)
lines(c(0, 100), c(1, 1), lty = 2, col = "red")

EnrichTests_RGB <- Compare2RangesWith3rd(L1CatalogGR, L1FragmGR, SegmGR, 
                                     SegmTable$itemRgb)
Colors <- sapply(names(EnrichTests_RGB$ORs), function(x) {
  eval(parse(text = paste("rgb(", x,", maxColorValue = 255)")))
})
barplot(EnrichTests_RGB$ORs, col = Colors)
lines(c(0, 100), c(1, 1), lty = 2, col = "red")

