# The script below analyzes the locations of LINE-1 insertions in the human
# genome, distingusihing between fragments, full-length and functional 
# insertions

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(ape)
library(seqinr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ShortRead)
library(csaw)

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file (Created in script AddColumns2L1Catalog.R)
L1CataloguePath <- "D:/L1polymORF/Data/L1CatalogExtended.csv"

# Path to chain file for b37 to hg19 (obtained from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/)
Chain_b37tohg19 <- "D:/L1polymORF/Data/b37tohg19.chain"

# Path to L1 GRanges from 1000 genome data
L1GRanges1000GenomesPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"

# Maximum L1 insertion length to be labeled a 'fragment'
MaxFragLength <- 5900

######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

# Read piRNA data
PiRNA_GR <- import.bed("D:/L1polymORF/Data/piR_hg19_sort.bed")

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable      <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")
RepeatTable_hg19 <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                     ranges = IRanges(start = RepeatTable$genoStart,
                                      end = RepeatTable$genoEnd),
                     strand = RepeatTable$strand)
L1RefGRFull  <- L1RefGR[width(L1RefGR) > 6000]
L1FragmGR    <- L1RefGR[width(L1RefGR) < MaxFragLength]
L1RefGR_hg19 <- GRanges(seqnames = RepeatTable_hg19$genoName,
                   ranges = IRanges(start = RepeatTable_hg19$genoStart,
                                    end = RepeatTable_hg19$genoEnd),
                   strand = RepeatTable_hg19$strand)
L1RefGRFull_hg19 <- L1RefGR_hg19[width(L1RefGR_hg19) > 6000]
L1FragmGR_hg19   <- L1RefGR_hg19[width(L1RefGR_hg19) < MaxFragLength]

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Lift catalog ranges to hg19
L1LiftoverList <- LiftoverL1Catalog(L1CatalogL1Mapped,  
                     ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR_hg19 <- L1LiftoverList$GRCatalogue_hg19

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
   ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                 L1CatalogL1Mapped$end_HG38),
                    end = pmax(L1CatalogL1Mapped$start_HG38,
                               L1CatalogL1Mapped$end_HG38)),
                      strand = L1CatalogL1Mapped$strand_L1toRef)

# Create genomic ranges for catalog L1 in the reference
blnInRef <- (L1CatalogL1Mapped$end_HG38 - L1CatalogL1Mapped$start_HG38) > 6000 
L1CatalogGR_Ref     <- L1CatalogGR[blnInRef]
L1CatalogGR_nonRef  <- L1CatalogGR[!blnInRef]
blnZero <- L1CatalogL1Mapped$Activity == 0
L1CatalogGR_nonZero <- L1CatalogGR[!blnZero]

# Get fullength L1 that are not in the catalog
blnOverlapCatalog <- overlapsAny(L1RefGRFull, L1CatalogGR, minoverlap = 6000)
L1RefGRFullnotCat <- L1RefGRFull[!blnOverlapCatalog]
L1GRZero          <- c(L1RefGRFullnotCat, L1CatalogGR[blnZero])

# Create genomic ranges for catalog L1 in the reference
blnInRef_hg19   <- width(L1CatalogGR_hg19) > 6000 
L1CatalogGR_Ref_hg19 <- L1CatalogGR_hg19[blnInRef_hg19]

# Get fullength L1 that are not in the catalog for hg19
blnOverlapCatalog_hg19 <- overlapsAny(L1RefGRFull_hg19, L1CatalogGR_hg19, 
                                      minoverlap = 6000)
L1RefGRFullnotCat_hg19 <- L1RefGRFull_hg19[!blnOverlapCatalog_hg19]

#######################################
#                                     #
#    Read & process 1000 genome data  #
#                                     #
#######################################

# Load data
load(L1GRanges1000GenomesPath)

# Subset genomic ranges to get full-length and fragments in three different 
# frequency classes
blnFull           <- L1_1000G_GR_hg19@elementMetadata@listData$InsLength >= 6000
GRL1Ins1000G_Full <- L1_1000G_GR_hg19[which(blnFull)]
FreqQuant_Full    <- quantile(GRL1Ins1000G_Full@elementMetadata@listData$Frequency)
blnAbove75        <- GRL1Ins1000G_Full@elementMetadata@listData$Frequency >=
                     FreqQuant_Full['75%']
GRL1Ins1000G_Full_HiF   <- GRL1Ins1000G_Full[blnAbove75]
GRL1Ins1000G_Full_lowF  <- GRL1Ins1000G_Full[!blnAbove75]
GRL1Ins1000G_Fragm      <- L1_1000G_GR_hg19[which(!blnFull)]

############################
#                          #
#    Basic statistics      #
#                          #
############################

##########
#  Count L1 per chromosome
##########

# L1 count per chromosome
L1CountPerChrom <- table(L1CatalogL1Mapped$Chromosome)
L1FragmPerChrom <- table(RepeatTable$genoName)

# Get chromosome length
load("D:/L1polymORF/Data/ChromLengthsHg38.RData")
ChromLengths <- ChromLengthsHg38
ChromNames   <- names(ChromLengthsHg38)

# Match chromosome names
ChromMatch          <- match(ChromNames, names(L1CountPerChrom))
ChromMatch_Fragm     <- match(ChromNames, names(L1FragmPerChrom))
L1CountMatched       <- L1CountPerChrom[ChromMatch]
L1CountMatched_Fragm <- L1FragmPerChrom[ChromMatch_Fragm]
names(L1CountMatched)[is.na(L1CountMatched)] <- ChromNames[is.na(L1CountMatched)]
names(L1CountMatched_Fragm)[is.na(L1CountMatched_Fragm)] <- 
  ChromNames[is.na(L1CountMatched_Fragm)]
L1CountMatched[is.na(L1CountMatched)] <- 0
L1CountMatched_Fragm[is.na(L1CountMatched_Fragm)] <- 0
  
# Fit poisson model
PoissonFit <- glm(L1CountMatched ~ ChromLengths, family = poisson(link = "identity"),
    start = c(0, 5*10^8))
summary(PoissonFit)
coefficients(PoissonFit)

# Plot chromosome length vs. L1 count
Cols <- rainbow(2)
plot(ChromLengths, L1CountMatched / max(L1CountMatched), 
     xlab = "Chromosome length", ylab = "L1 count", col = Cols[1], pch = 16)
points(ChromLengths, L1CountMatched_Fragm / max(L1CountMatched_Fragm), col = Cols[2], pch = 16)
legend("bottomright", col = Cols, legend = c("Full-length", "Fragments"), pch = c(16,16))
#CreateDisplayPdf('D:/L1polymORF/Figures/L1NrInsertsVsChromLength.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1NrInsertsVsChromLength.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Plot chromosome length vs. L1 fragment count
plot(ChromLengths, L1CountMatched_Fragm, xlab = "Chromosome length", 
     ylab = "L1 fragment count")

##########
#  Count L1 intersecting with exons, genes and promoters
##########

# Get ranges of genes
GRgenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
GRgenes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Count overlaps between L1 and exons 
blnL1OverlapExon     <- overlapsAny(L1CatalogGR, 
                               exons(TxDb.Hsapiens.UCSC.hg38.knownGene)) 
blnExonOverlapL1     <- overlapsAny(exons(TxDb.Hsapiens.UCSC.hg38.knownGene), L1CatalogGR) 
sum(blnL1OverlapExon)

# Count overlaps between L1 and exons 
blnL1OverlapGene_Cat   <- overlapsAny(L1CatalogGR, GRgenes)
blnL1OverlapGene_Fragm <- overlapsAny(L1FragmGR, GRgenes)
blnGeneOverlapCat    <- overlapsAny(GRgenes, L1CatalogGR)
blnGeneOverlapFragm  <- overlapsAny(GRgenes, L1FragmGR)
blnL1OverlapPromoter_Cat <- overlapsAny(L1CatalogGR, 
   promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 5000)) 
mean(blnL1OverlapGene_Cat)
mean(blnL1OverlapGene_Fragm)
sum(blnL1OverlapPromoter_Cat)
mean(blnL1OverlapGene_Cat)

# Perform logistic regression to determine whether fragment length predicts
# probability of gene overlap
blnL1OverlapGene_Fragm <- overlapsAny(L1FragmGR, GRgenes)
FitData <- data.frame(blnOverlap = blnL1OverlapGene_Fragm,
                      FragmWidth = width(L1FragmGR))
LogRegFit <- glm(blnOverlap ~ FragmWidth, data = FitData,
                 family = binomial)
summary(LogRegFit)
par(mfrow = c(1, 1))
plot(c(0, 3), c(0, 0.2), type = "n", xlab = "L1 type", ylab = "Proportion in gene",
     xaxt = "n")
axis(1, at = 1:2, labels = c("Fragment", "Full-length"))
MeanProps <- c(mean(blnL1OverlapGene_Fragm), mean(blnL1OverlapGene_Cat))
StErr <- sqrt(c(var(blnL1OverlapGene_Fragm), var(blnL1OverlapGene_Cat))/
                    c(length(L1FragmGR), length(L1CatalogGR)))
rect(c(0.7, 1.7), c(0, 0), c(1.3, 2.3), MeanProps, col = "grey")
AddErrorBars(MidX = c(1, 2), MidY = MeanProps, ErrorRange = StErr,
             TipWidth = 0.1)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntron.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntron.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Calculate proportions of L1 insertions intersecting with gene bodies
mean(blnL1OverlapGene_Cat)
mean(blnL1OverlapGene_Fragm)
pbinom(q = sum(blnL1OverlapGene_Cat), size  = length(blnL1OverlapGene_Cat),
       prob = mean(blnL1OverlapGene_Fragm))

# Plot proportion overlap per insert size class 
WidthCut <- cut(width(L1FragmGR), breaks = seq(0, 5000, 500))
PropInGene <- aggregate(blnL1OverlapGene_Fragm ~ WidthCut, FUN = mean)
VarInGene <- aggregate(blnL1OverlapGene_Fragm ~ WidthCut, FUN = var)
NInGene <- aggregate(blnL1OverlapGene_Fragm ~ WidthCut, FUN = length)
StErrInGene <- sqrt(VarInGene$blnL1OverlapGene_Fragm / 
                      NInGene$blnL1OverlapGene_Fragm)
AvWidth    <- aggregate(width(L1FragmGR) ~ WidthCut, FUN = mean)
plot(AvWidth$`width(L1FragmGR)`, PropInGene$blnL1OverlapGene_Fragm,
     xlab = "Fragment size [bp]", ylab = "Proportion in introns",
     pch = 16, type = "n", ylim = c(0.05, 0.25))
AddBars(Ys = PropInGene$blnL1OverlapGene_Fragm, MidX = AvWidth$`width(L1FragmGR)`)
  
AddErrorBars(MidX = AvWidth$`width(L1FragmGR)`, MidY = PropInGene$blnL1OverlapGene_Fragm, 
             ErrorRange = StErrInGene,TipWidth = 100)
WidthOrder <- order(width(L1FragmGR))
lines(width(L1FragmGR)[WidthOrder], fitted.values(LogRegFit)[WidthOrder])
segments(0, mean(blnL1OverlapGene_Cat), 10^6, mean(blnL1OverlapGene_Cat), col = "red",
         lty = 2)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntronVsSize.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntronVsSize.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Get predicted logit transformed proportion and standard error for prediction
# for length of 6000 bp
NewData  <- data.frame(FragmWidth = 6064)
PredictP <- predict(LogRegFit, se.fit = T, newdata = NewData)

# Sample proportion and number of intersecting L1
SampleSize  <- 10000
SampleFit   <- rnorm(SampleSize, PredictP$fit, PredictP$se.fit)
SampleProps <- exp(SampleFit) / (1 + exp(SampleFit))
hist(SampleProps)
SampledVals <- sapply(SampleProps, function(x) {
  rbinom(1, length(L1CatalogGR), x)
})
ObsVal <- sum(blnL1OverlapGene_Cat)
sum(SampledVals >= ObsVal) / length(SampledVals)
hist(SampledVals, main = "", xlab = "Number of L1 in gene bodies")
segments(x0 = ObsVal, y0 = 0, y1 = 5000, col = "red")
#CreateDisplayPdf('D:/L1polymORF/Figures/L1propInL1.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1propInL1.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

#1 - pbinom(sum(blnL1OverlapGene), length(L1CatalogGR), PredictProp)

##################################
#                                #
#    Distance distributions      #
#                                #
##################################

##########
#  Defining auxiliary function
##########

# Auxiliary function to get distances to closest gene
Dist2ClosestGene <- function(GR){
  DistGeneObj <- distanceToNearest(GR, GRgenes, ignore.strand = T) 
  DistGeneObj@elementMetadata@listData$distance
}
Dist2Closest <- function(GR1, GR2){
  DistObj <- distanceToNearest(GR1, GR2, ignore.strand = T) 
  DistObj@elementMetadata@listData$distance
}

# Auxiliary unction to create genomic ranges for both sides of a loop
getLoopGRs <- function(FileName){
  Loops <- read.delim(FileName)
  Loops$chr1 <- paste("chr", Loops$chr1, sep = "")
  Loops$chr2 <- paste("chr", Loops$chr2, sep = "")
  LoopsGR1 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr1", 
                                       start.field="x1", end.field = "x2")
  LoopsGR2 <- makeGRangesFromDataFrame(Loops, seqnames.field = "chr2", 
                                       start.field="y1", end.field = "y2")
  blnOverlapLoop <- overlapsAny(LoopsGR1, LoopsGR2)
  AllLoops <- c(LoopsGR1, LoopsGR2[!blnOverlapLoop])
  GRangesList(LoopsGR1 = LoopsGR1, LoopsGR2 = LoopsGR2, AllLoops = AllLoops)
}

# Calculate distances from fragment L1, catalog L1 in the reference genome and
# full-length L1 that are not in the catalog
L1DistGene_Fragm      <- Dist2ClosestGene(L1FragmGR)
L1DistGene_CatRef     <- Dist2ClosestGene(L1CatalogGR_Ref)
L1DistGene_FullnotCat <- Dist2ClosestGene(L1RefGRFullnotCat)

# Calculate distances from fragment L1, catalog L1 in the reference genome and
# full-length L1 that are not in the catalog
L1DistGene_Fragm_hg19      <- Dist2Closest(L1FragmGR_hg19, GRgenes_hg19)
L1DistGene_CatRef_hg19     <- Dist2Closest(L1CatalogGR_Ref_hg19, GRgenes_hg19)
L1DistGene_FullnotCat_hg19 <- Dist2Closest(L1RefGRFullnotCat_hg19, GRgenes_hg19)

# Calculate distances from fragment L1, catalog L1 in the reference genome and
# full-length L1 that are not in the catalog to piRNA
L1DistPiRNA_Fragm_hg19      <- Dist2Closest(L1FragmGR_hg19, PiRNA_GR)
L1DistPiRNA_CatRef_hg19     <- Dist2Closest(L1CatalogGR_Ref_hg19, PiRNA_GR)
L1DistPiRNA_FullnotCat_hg19 <- Dist2Closest(L1RefGRFullnotCat_hg19, PiRNA_GR)

# Calculate distances from full-length and fragment L1 from 1000 genome data to
# closest gene
L1DistGene_1000GFull      <- Dist2Closest(GRL1Ins1000G_Full, GRgenes_hg19)
L1DistGene_1000GFullLow   <- Dist2Closest(GRL1Ins1000G_Full_lowF, GRgenes_hg19)
L1DistGene_1000GFullHigh  <- Dist2Closest(GRL1Ins1000G_Full_HiF, GRgenes_hg19)
L1DistGene_1000GFragm     <- Dist2Closest(GRL1Ins1000G_Fragm, GRgenes_hg19)

##########
#  Calculate distances to domain boundaries, catalog L1
##########

# List of files with domains and loops
DomainFiles <- list.files("D:/L1polymORF/Data/HiCData/", pattern = "domainlist.txt",
                          full.names = T)
DomainFiles <- DomainFiles[! DomainFiles %in% grep(".gz", DomainFiles, value = T)]
LoopFiles <- list.files("D:/L1polymORF/Data/HiCData/", pattern = "looplist.txt",
                          full.names = T)
LoopFiles <- LoopFiles[! LoopFiles %in% grep(".gz", LoopFiles, value = T)]

# Get list of domain and loop ranges
DomainGRList   <- lapply(DomainFiles, getLoopGRs)
LoopGRList     <- lapply(LoopFiles, getLoopGRs)
names(DomainGRList) <- sapply(DomainFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
names(LoopGRList) <- sapply(LoopFiles, 
                              function(x) strsplit(x, "_")[[1]][2])
DomainDistList <- lapply(DomainGRList, function(x) {
  list(DistFragm = Dist2Closest(L1FragmGR_hg19, x$AllLoops),
       DistCatRef = Dist2Closest(L1CatalogGR_Ref_hg19, x$AllLoops),
       DistRefnotCat = Dist2Closest(L1RefGRFullnotCat_hg19, x$AllLoops))
})
DomainDistList_1000G <- lapply(DomainGRList, function(x) {
  list(DistFragm = Dist2Closest(GRL1Ins1000G_Fragm, x$AllLoops),
       DistFullHigh = Dist2Closest(GRL1Ins1000G_Full_HiF, x$AllLoops),
       DistFullLow = Dist2Closest(GRL1Ins1000G_Full_lowF, x$AllLoops))
})

LoopDistList <- lapply(LoopGRList, function(x) {
  list(DistFragm = Dist2Closest(L1FragmGR_hg19, x$AllLoops),
       DistCatRef = Dist2Closest(L1CatalogGR_Ref_hg19, x$AllLoops),
       DistRefnotCat = Dist2Closest(L1RefGRFullnotCat_hg19, x$AllLoops))
})
DomainDist0 <- sapply(DomainDistList, function(x) {
  c(nrFragm = length(x$DistFragm),
    nrCatRef = length(x$DistCatRef),
    nrRefnotCat = length(x$DistRefnotCat),
    prop0Fragm = mean(x$DistFragm == 0),
    prop0CatRef = mean(x$DistCatRef == 0),
    prop0RefnotCat = mean(x$DistRefnotCat == 0),
    sum0Fragm = sum(x$DistFragm == 0),
    sum0CatRef = sum(x$DistCatRef == 0),
    sum0RefnotCat = sum(x$DistRefnotCat == 0))
})
LoopDist0 <- sapply(LoopDistList, function(x) {
  c(nrFragm = length(x$DistFragm),
    nrCatRef = length(x$DistCatRef),
    nrRefnotCat = length(x$DistRefnotCat),
    prop0Fragm = mean(x$DistFragm == 0),
    prop0CatRef = mean(x$DistCatRef == 0),
    prop0RefnotCat = mean(x$DistRefnotCat == 0),
    sum0Fragm = sum(x$DistFragm == 0),
    sum0CatRef = sum(x$DistCatRef == 0),
    sum0RefnotCat = sum(x$DistRefnotCat == 0))
})

# Probability of observed number of catalog intersections given the proportion
# of fragment interesections
apply(LoopDist0, 2, FUN = function(x){
  1 - pbinom(x["sum0CatRef"] - 1, x["nrCatRef"], x["prop0Fragm"])
})
apply(DomainDist0, 2, FUN = function(x){
  1 - pbinom(x["sum0CatRef"] - 1, x["nrCatRef"], x["prop0Fragm"])
})


# # Create different genomic ranges
# LoopsGR_GM12878Dom <- getLoopGRs("D:/L1polymORF/Data/HiCData/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt")
# Loops_HMECDom      <- getLoopGRs("D:/L1polymORF/Data/HiCData/GSE63525_HMEC_Arrowhead_domainlist.txt")
# Loops_GM12878Loop  <- getLoopGRs("D:/L1polymORF/Data/HiCData/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt")
# 
# DistFragm <- Dist2Closest(L1FragmGR_hg19, AllLoops)
# DistCat   <- Dist2Closest(L1CatalogGR_hg19, AllLoops)
# DistFragm <- Dist2Closest(GRL1Ins1000G_Fragm, AllLoops)
# DistCat   <- Dist2Closest(GRL1Ins1000G_Full_HiF, AllLoops)
# DistCatLF   <- Dist2Closest(GRL1Ins1000G_Full_lowF, AllLoops)

############################
#                          #
#    QQ plots              #
#                          #
############################

##########
#  Auxilliary function to create qq plot to compare distance distributions
##########

QQDistPlot <- function(FragmDist, Dist1, Dist2 = NULL, Dist3 = NULL, 
                       xLab = "", 
                       yLab = "", 
                       NrSamples = 10000,
                       QuantV = seq(0, 1, 0.01),
                       Title = "", blnAxLab = T,
                       XLim = NULL, YLim = NULL){
  QSampled <- SampleQuantiles(FragmDist, length(Dist1),
                              QuantV = QuantV, NrSamples = NrSamples)
  QSMat <- QSampled$QMat
  idxF <- 1:ncol(QSMat)
  idxR <- ncol(QSMat):1
  QQ1  <- qqplot(FragmDist, Dist1, plot.it = F)
  AllDist <- c(Dist1, Dist2, Dist3)
  if (is.null(XLim)) XLim <-  c(min(FragmDist), max(FragmDist))
  if (is.null(YLim)) YLim <-  c(min(AllDist), max(AllDist))
  plot(QQ1$x, QQ1$y, xlab = xLab, ylab = yLab, main = Title, 
       xlim = XLim, ylim = YLim)
  polygon(QSMat[1,c(idxF, idxR)], c(QSMat[2, ], QSMat[3, idxR]), 
          col = "grey", border = NA)
  points(QQ1$x, QQ1$y)
  if (!is.null(Dist2)){
    QQ2 <- qqplot(FragmDist, Dist2, plot.it = F)
    points(QQ2$x, QQ2$y, pch = 2)
  }
  if (!is.null(Dist3)){
    QQ3 <- qqplot(FragmDist, Dist3, plot.it = F)
    points(QQ3$x, QQ3$y, pch = 3)
  }
  lines(c(0, 10^10), c(0, 10^10))
  mean(QSampled$SampleMeans <= mean(Dist1))
}

# QQplot for different cell lines
par(mfrow = c(3, 3), mar = c(3, 2, 2, 0.5), oma = c(2, 2.5, 1, 1))
QSamplesDomain <- sapply (1:length(DomainDistList), function(i){
  x <- DomainDistList[[i]] 
  QQDistPlot(x$DistFragm, x$DistCatRef, x$DistRefnotCat, 
             Title = names(DomainDistList)[i])
})
names(QSamplesDomain) <- names(DomainDistList)
mtext("Distance fragment L1 to domain", side = 1, line = 1, outer = T)
mtext("Distance full-length L1 to domain", side = 2, line = 1, outer = T)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1DomainDistQQ.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1DomainDistQQ.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

par(mfrow = c(3, 3), mar = c(3, 2, 2, 0.5), oma = c(2, 2.5, 1, 1))
QSamplesDomain_1000G <- sapply (1:length(DomainDistList_1000G), function(i){
  x <- DomainDistList_1000G[[i]] 
  QQDistPlot(x$DistFragm, x$DistFullHigh, x$DistFullLow, 
             Title = names(DomainDistList_1000G)[i])
})
names(QSamplesDomain_1000G) <- names(DomainDistList_1000G)
mtext("Distance fragment L1 to domain", side = 1, line = 1, outer = T)
mtext("Distance full-length L1 to domain", side = 2, line = 1, outer = T)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1DomainDistQQ_1000G.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1DomainDistQQ_1000G.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

plot(QSamplesDomain, QSamplesDomain_1000G)

par(mfrow = c(3, 3)) 
QSamplesLoop <- sapply (1:length(LoopDistList), function(i){
    x <- LoopDistList[[i]] 
  QQDistPlot(x$DistFragm, x$DistCatRef, x$DistRefnotCat, 
             Title = names(DomainDistList)[i])
})
names(QSamplesLoop) <- names(LoopDistList)
mtext("Distance fragment L1 to loop borders", side = 1, line = 1, outer = T)
mtext("Distance full-length L1 to loop borders", side = 2, line = 1, outer = T)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1LoopDistQQ.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1LoopDistQQ.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

par(mfrow = c(1, 1))
plot(QSamplesLoop, QSamplesDomain)

# Plot quantiles against each other for distance to genes and distance to domains
par(mfrow = c(2, 2), mar = c(3, 2, 2, 0.5), oma = c(2, 2.5, 1, 1))
YI <- 0.2
pGeneDistCat <- QQDistPlot(L1DistGene_Fragm, L1DistGene_CatRef, 
                           L1DistGene_FullnotCat, Title = "Distance to genes")
legend("bottomright", legend = c("intact", "not intact"), pch = c(1,2),
       y.intersp = YI, bty = "n")
# pGeneDist1000G <- QQDistPlot(L1DistGene_1000GFragm, L1DistGene_1000GFullHigh, 
#                              L1DistGene_1000GFullLow)
x <- DomainDistList[[1]] 
pDomainDistCat <- QQDistPlot(x$DistFragm, x$DistCatRef, x$DistRefnotCat,
                             Title = "Distance to loops GM12878")
pDomainDistCat <- QQDistPlot(x$DistFragm, x$DistCatRef, x$DistRefnotCat,
                             Title = "Distance to loops GM12878",
                             XLim = c(0, 2*10^6), YLim = c(0, 2*10^6))
pDomainDistCat <- QQDistPlot(x$DistFragm, x$DistCatRef, x$DistRefnotCat,
                             Title = "Distance to loops GM12878",
                             XLim = c(0, 10^5), YLim = c(0, 10^5))
mtext("Distance from fragment L1", side = 1, line = 1, outer = T)
mtext("Distance full-length L1", side = 2, line = 1, outer = T)
#CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQQ_Catalog.pdf')
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQQ_Catalog.pdf', 
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"')

# Plot quantiles against each other for distance to piRNA and distance to domains
par(mfrow = c(1,1), mar = c(3, 2, 2, 0.5), oma = c(2, 2.5, 1, 1))
YI <- 0.2
piRNADistCat <- QQDistPlot(L1DistPiRNA_Fragm_hg19, L1DistPiRNA_CatRef_hg19, 
                           Title = "Distance to piRNA")
legend("bottomright", legend = c("intact", "not intact"), pch = c(1,2),
       y.intersp = YI, bty = "n")


##################################
#                                #
#    Other tests      #
#                                #
##################################

# Get index of nearest gene
idxNearest <- nearest(L1CatalogGR, GRgenes, ignore.strand = T)

# Test whether there is a relationship between strandedness of L1 and the
# nearest gene
NearestStrand <- as.vector(strand(GRgenes)[idxNearest])
table(NearestStrand, as.vector(strand(L1CatalogGR)))
chisq.test(NearestStrand, as.vector(strand(L1CatalogGR)))

# Test whether there is a relationship between strandedness of L1 and the
# nearest gene
UpStream   <- start(GRgenes)[idxNearest] < start(L1CatalogGR) 
table(UpStream, as.vector(strand(L1CatalogGR)))
chisq.test(UpStream, as.vector(strand(L1CatalogGR)))

# Perform regression to determine whether distance to closest gene depends on
# fragment length
FitDistVsLength <- glm(L1DistGene_Fragm ~ width(L1FragmGR))
summary(FitDistVsLength)

mean(L1DistPiRNA_Fragm_hg19 == 0)
Nr
pbinom(length(L1DistPiRNA_CatRef_hg19))

QQDistPlot(L1DistGene_1000GFragm, L1DistGene_1000GFullHigh, L1DistGene_1000GFullLow)
QQDistPlot(L1DistGene_1000GFullLow, L1DistGene_1000GFullHigh)
QQDistPlot(L1DistGene_1000GFragm, L1DistGene_1000GFull, 
           xLab = "Distance from L1 fragments",
           yLab = "Distance from full-length L1")

