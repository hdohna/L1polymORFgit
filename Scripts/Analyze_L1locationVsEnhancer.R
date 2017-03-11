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
RegTable      <- read.table("D:/L1polymORF/Data/ChromHMM", header = T)
EnhancerTable <- RegTable[grep("Enhancer", RegTable$name), ]
unique(EnhancerTable$name)
# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
RegGR <- makeGRangesFromDataFrame(RegTable, start.field = "ChromStart", 
                                       end.field = "ChromEnd")
idxSmall <- width(RegGR) <= 10000
RegGR <- RegGR[idxSmall]
RegTable <- RegTable[idxSmall,]
EnhancerGR <- makeGRangesFromDataFrame(EnhancerTable, start.field = "ChromStart", 
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
#    Calculate distances to enhancers and genes  #
#                                                #
##################################################

# Auxiliary function to get distances to closest gene
Dist2ClosestGR <- function(GR1, GR2){
  DistObj <- distanceToNearest(GR1, GR2, ignore.strand = T) 
  DistObj@elementMetadata@listData$distance
}

# Get ranges of genes
GRgenes         <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
GenesGR_reduced <- GRanges(seqnames = seqnames(GRgenes), ranges(GRgenes))
GRcombined      <- c(GenesGR_reduced, EnhancerGR)

# Calculate distances from full-length L1 to nearest feature
L1Dist_cat <- Dist2ClosestGR(L1CatalogGR, EnhancerGR)
L1Dist_cat_Combined <- Dist2ClosestGR(L1CatalogGR, GRcombined)
hist(L1CatDist)

idxNearestCat     <- nearest(L1CatalogGR, RegGR, ignore.strand = T)
OverlapTableCat   <- RegTable[idxNearestCat, ]
idxNearestFragm   <- nearest(L1FragmGR, RegGR, ignore.strand = T)
OverlapTableFragm <- RegTable[idxNearestFragm, ]
FreqCat   <- table(OverlapTableCat$name)
FreqFragm <- table(OverlapTableFragm$name)
chisq.test(FreqCat, FreqFragm)
plot(as.numeric(FreqCat), as.numeric(FreqFragm))
text(as.numeric(FreqCat), as.numeric(FreqFragm),names(FreqCat))
lines(c(0, 1000) * length(L1CatalogGR), c(0, 1000) * length(L1FragmGR))

# Calculate distances from fragment L1 to nearest enhancer
idxNearestCat     <- nearest(L1CatalogGR, EnhancerGR, ignore.strand = T)
OverlapTableCat   <- EnhancerTable[idxNearestCat, ]
idxNearestFragm   <- nearest(L1FragmGR, EnhancerGR, ignore.strand = T)
OverlapTableFragm <- EnhancerTable[idxNearestFragm, ]
FreqCat   <- table(OverlapTableCat$name)
FreqFragm <- table(OverlapTableFragm$name)
chisq.test(FreqCat, FreqFragm)
plot(as.numeric(FreqCat), as.numeric(FreqFragm))
text(as.numeric(FreqCat), as.numeric(FreqFragm),names(FreqCat))
lines(c(0, 1000) * length(L1CatalogGR), c(0, 1000) * length(L1FragmGR))



L1Dist_Fragm <- Dist2ClosestGR(L1FragmGR, EnhancerGR)
Poverlap <- sum(L1Dist_Fragm == 0) / length(L1Dist_Fragm)
sum(L1Dist_cat == 0) / length(L1Dist_cat)
Poverlap^length(L1Dist_cat)
L1Dist_Fragm <- Dist2ClosestGR(L1FragmGR, GRcombined)

# Get the distance closest gene for full-length L1 that are not in the catalog 
L1DistGene_FullnotCat <- Dist2ClosestGene(L1RefGRFullnotCat)


# Plot quantiles against each other
par(mfrow = c(1, 1))
qqplot(L1Dist_Fragm, L1Dist_cat, ylab = "Distance full-length L1 to gene",
       xlab = "Distance fragment L1 to gene")
lines(c(0, 10^7), c(0, 10^7))
QQ2 <- qqplot(L1DistGene_Fragm, L1RefDistGene, plot.it = F)
points(QQ2$x, QQ2$y, pch = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQQ_Catalog.pdf')

# Plot quantiles against each other comparing catalog and non-catalog full
# length L1
L1NonZeroDistGene <- Dist2ClosestGene(L1CatalogGR_nonZero)
L1ZeroDistGene    <- Dist2ClosestGene(L1GRZero)
par(mfrow = c(1, 1))
#qqplot(L1ZeroDistGene, L1NonZeroDistGene,
L1CatalogL1Mapped$Accession[which(L1DistGene > 1500000)]
qqplot(L1DistGene_FullnotCat, L1DistGene,
              ylab = "Distance catalog full-length L1 to gene",
       xlab = "Distance non-catalog full-length L1 to gene")
lines(c(0, 10^7), c(0, 10^7))
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQQ_CatalogVsNonCatalog.pdf')


# Calculate distances from full-length and fragment L1 from 1000 genome data to
# closest gene
L1Full1000G_DistGene    <- Dist2ClosestGene(GRL1Ins1000G_Full)
L1FullMed1000G_DistGene <- Dist2ClosestGene(GRL1Ins1000G_Full_MedF)
L1FullHi1000G_DistGene  <- Dist2ClosestGene(GRL1Ins1000G_Full_HiF)

L1Fragm1000G_DistGene    <- Dist2ClosestGene(GRL1Ins1000G_Fragm)
L1FragmMed1000G_DistGene <- Dist2ClosestGene(GRL1Ins1000G_Fragm_MedF)
L1FragmHi1000G_DistGene  <- Dist2ClosestGene(GRL1Ins1000G_Fragm_HiF)

qqplot(L1Fragm1000G_DistGene, L1Full1000G_DistGene, ylab = "Distance full-length L1 to gene",
       xlab = "Distance fragment L1 to gene")
lines(c(0, 10^7), c(0, 10^7))
# QQMed <- qqplot(L1FragmMed1000G_DistGene, L1FullMed1000G_DistGene, plot.it = F)
# points(QQMed$x, QQMed$y, pch = 2)
QQHi <- qqplot(L1FragmHi1000G_DistGene, L1FullHi1000G_DistGene, plot.it = F)
points(QQHi$x, QQHi$y, pch = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQQ_1000G.pdf')



# Plot histograms for catalog and fragment distance distribution
Hist_Fragm <- hist(L1DistGene_Fragm, breaks = seq(0, 5*10^6, 5*10^4),
                   plot = F)
Hist_Catalog <- hist(L1DistGene, breaks = seq(0, 5*10^6, 5*10^4),
                   plot = F)

par(mfrow = c(1, 1))
plot(Hist_Catalog$mids, log(Hist_Catalog$density), ylab = "Density", 
     xlab = "Distance to closest gene", col = "blue", ylim = c(-20, -10))
points(Hist_Fragm$mids, log(Hist_Fragm$density), col = "red")
legend("topright", legend = c("catalog", "fragment"), col = c("blue", "red"),
       lty = c(1,1))
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistLogHisto.pdf')

##########
#  Compare distances to genes with other properties
##########

plot(L1DistGene, L1CatalogL1Mapped$Activity)
plot(L1DistGene, L1CatalogL1Mapped$Allele_frequency)

