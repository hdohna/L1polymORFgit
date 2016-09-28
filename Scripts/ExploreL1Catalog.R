# The script below explore a catalog of full-length L1

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
library(ShortRead)
library(csaw)
library(quantreg)
library(car)
# library(GO.db)
# library(org.Hs.eg.db)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
#L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_16-33-31_2016.csv"

# Path to L1 GRanges from 1000 genome data
L1GRanges1000GenomesPath <- "D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData"

# Path to bam file of mapped L1 fragments
BamFilePath <- 'D:/L1polymORF/Data/L1fragments_aln2consens.sorted.bam'

# Path to L1 alignment file 
AlignFileName <- 'D:/L1polymORF/Data/L1HSSequences_L100_withConsens_aligned.fas'

# Length of Flanking sequence to be used for alignment
FlankLength <- 200

# Maximum fragment length
MaxFragLength <- 5900


######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

# Get reads from bam file aligned to L1HS concensus
L1ConsensGR <- GRanges(seqnames = 'L1HS_L1_Homo_sapiens',
                       IRanges(start = 1, end = 6064))
readParams  <- ScanBamParam(which = L1ConsensGR,
                            what = c('qname', "pos","cigar", "qwidth", "seq"), tag = "NM")
Reads       <- scanBam(BamFilePath, param = readParams)
all(Reads[[1]]$qwidth == nchar(Reads[[1]]$seq))
hist(Reads[[1]]$qwidth)
Dist2Consens <- Reads[[1]]$tag$NM / Reads[[1]]$qwidth
names(Dist2Consens) <- Reads[[1]]$qname
hist(Dist2Consens)
Reads[[1]]$seq[1:10]
Reads[[1]]$qname[1:10]

# Match pattern following ACA motif at 3' end of L1
MatchList <- vmatchPattern('TTAGTGGGTG', Reads[[1]]$seq)
NrMatches <- sapply(MatchList, length)
table(NrMatches)
ACAmotif <- sapply(1:length(MatchList), function(i){
  if (length(MatchList[[i]]) > 0){
    End <- max(1, start(MatchList[[i]])[1] - 1)
    as.character(subseq(Reads[[1]]$seq[i], max(1, End - 2), End))
  } else {
    NA
  }
})
table(ACAmotif)
names(ACAmotif) <- Reads[[1]]$qname

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
NameMatch <- match(L1FragmNames, names(Dist2Consens))
Dist2Consens <- Dist2Consens[NameMatch]
hist(Dist2Consens)
cor(width(L1FragmGR), Dist2Consens)
sum(Dist2Consens < 0.4, na.rm = T)

# Get 2 sets of fragments: close to the L1HS consensus and divergent
L1FragmGR_close   <- L1FragmGR[which(Dist2Consens < 0.02)]
L1FragmGR_diverge <- L1FragmGR[which(Dist2Consens >= 0.02)]

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)
L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create an overview table of L1 counts
# table(L1Catalogue$Reference)
# table(blnInRef1[blnAllele1], L1Catalogue$Reference[blnAllele1])
# L1Catalogue[L1Catalogue$Reference == "Brouha2003" & !blnInRef1, 1:5]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
   ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                 L1CatalogL1Mapped$end_HG38),
                    end = pmax(L1CatalogL1Mapped$start_HG38,
                               L1CatalogL1Mapped$end_HG38)),
                      strand = L1CatalogL1Mapped$strand_L1toRef)

# Create genomic ranges for catalog L1 in the reference
blnInRef <- (L1CatalogL1Mapped$end_HG38 - L1CatalogL1Mapped$start_HG38) > 6000 
L1CatalogGR_Ref <- L1CatalogGR[blnInRef]
L1CatalogGR_nonRef <- L1CatalogGR[!blnInRef]
blnZero <- L1CatalogL1Mapped$Activity == 0
L1CatalogGR_nonZero <- L1CatalogGR[!blnZero]


# Get fullength L1 that are not in the catalog
blnOverlapCatalog <- overlapsAny(L1RefGRFull, L1CatalogGR, minoverlap = 6000)
L1RefGRFullnotCat <- L1RefGRFull[!blnOverlapCatalog]
L1GRZero <- c(L1RefGRFullnotCat, L1CatalogGR[blnZero])

##########
#  Load genomic ranges from 1000 genome data
##########

# Load data
load(L1GRanges1000GenomesPath)

# Subset L1Ins1000G_Info so that it only contains the infor for L1 mapped to hg38
L1Ins1000G_Info_mapped <- L1Ins1000G_Info[idxUniqueMapped_hg38,]

# Subset genomic ranges to get full-length and fragments in three different 
# frequency classes
blnFull    <- L1Ins1000G_Info_mapped$InsLength >= 6000
blnMedFreq <- L1Ins1000G_Info_mapped$Freq >= 10
blnHiFreq  <- L1Ins1000G_Info_mapped$Freq >= 20
GRL1Ins1000G_Full       <- GRL1Ins1000G_hg38Mapped[which(blnFull)]
GRL1Ins1000G_Full_MedF  <- GRL1Ins1000G_hg38Mapped[which(blnFull & blnMedFreq)]
GRL1Ins1000G_Full_HiF   <- GRL1Ins1000G_hg38Mapped[which(blnFull & blnHiFreq)]
GRL1Ins1000G_Fragm      <- GRL1Ins1000G_hg38Mapped[which(!blnFull)]
GRL1Ins1000G_Fragm_MedF <- GRL1Ins1000G_hg38Mapped[which(!blnFull & blnMedFreq)]
GRL1Ins1000G_Fragm_HiF  <- GRL1Ins1000G_hg38Mapped[which(!blnFull & blnHiFreq)]

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
ChromLengths <-ChromLengthsHg38

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
CreateDisplayPdf('D:/L1polymORF/Figures/L1NrInsertsVsChromLength.pdf')

# Plot chromosome length vs. L1 fragment count
plot(ChromLengths, L1CountMatched_Fragm, xlab = "Chromosome length", 
     ylab = "L1 fragment count")

##########
#  Count L1 intersecting with exons, genes and promoters
##########

# Get ranges of genes
GRgenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Count overlaps with exons and promoters
blnL1OverlapExon     <- overlapsAny(L1CatalogGR, 
                               exons(TxDb.Hsapiens.UCSC.hg38.knownGene)) 
L1CatalogGR[blnL1OverlapExon]
blnExonOverlapL1     <- overlapsAny(exons(TxDb.Hsapiens.UCSC.hg38.knownGene), L1CatalogGR) 
L1CatalogGR[blnL1OverlapExon]
blnL1OverlapGene     <- overlapsAny(L1CatalogGR, GRgenes)
blnGeneOverlapCat    <- overlapsAny(GRgenes, L1CatalogGR)
blnGeneOverlapFragm  <- overlapsAny(GRgenes, L1FragmGR)
blnL1OverlapPromoter <- overlapsAny(L1CatalogGR, 
   promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 5000)) 
mean(blnGeneOverlapCat)
mean(blnGeneOverlapCat)
sum(blnL1OverlapExon)
sum(blnL1OverlapGene)
sum(blnL1OverlapPromoter)
mean(blnL1OverlapGene)

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
MeanProps <- c(mean(blnL1OverlapGene_Fragm), mean(blnL1OverlapGene))
StErr <- sqrt(c(var(blnL1OverlapGene_Fragm), var(blnL1OverlapGene))/
                    c(length(L1FragmGR), length(L1CatalogGR)))
rect(c(0.7, 1.7), c(0, 0), c(1.3, 2.3), MeanProps, col = "grey")
AddErrorBars(MidX = c(1, 2), MidY = MeanProps, ErrorRange = StErr,
             TipWidth = 0.1)
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntron.pdf')

# Calculate proportions of L1 insertions intersecting with gene bodies
mean(blnL1OverlapGene)
mean(blnL1OverlapGene_Fragm)
pbinom(q = sum(blnL1OverlapGene), size  = length(blnL1OverlapGene),
       prob = mean(blnL1OverlapGene_Fragm))

# Get proportion overlap 
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
segments(0, mean(blnL1OverlapGene), 10^6, mean(blnL1OverlapGene), col = "red",
         lty = 2)
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntronVsSize.pdf')

# Test whether the observed proportion of catalog L1 overlapping with genes is
# statistically significant from the proportion expected at size 5000
idxMax <- which.max(width(L1FragmGR))
PredictProp <- fitted.values(LogRegFit)[idxMax]

# Calculate standard error for prediction
NewData <- data.frame(FragmWidth = 6000)

PredictP <- predict(LogRegFit, se.fit = T, newdata = NewData)
PredictP$se.fit[idxMax]

# Sample proportion and number of intersecting L1
SampleSize  <- 10000
SampleFit <- rnorm(SampleSize, PredictP$fit, PredictP$se.fit)
SampleProps <- exp(SampleFit) / (1 + exp(SampleFit))
hist(SampleProps)
SampledVals <- sapply(SampleProps, function(x) {
  rbinom(1, length(L1CatalogGR), x)
})
ObsVal <- sum(blnL1OverlapGene)
sum(SampledVals >= ObsVal) / length(SampledVals)
hist(SampledVals, main = "", xlab = "Number of L1 in gene bodies")
segments(x0 = ObsVal, y0 = 0, y1 = 5000, col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1propInL1.pdf')

#1 - pbinom(sum(blnL1OverlapGene), length(L1CatalogGR), PredictProp)

##########
#  Calculate distances to genes
##########

# Auxiliary function to get distances to closest gene
Dist2ClosestGene <- function(GR){
  DistGeneObj <- distanceToNearest(GR, GRgenes, ignore.strand = T) 
  DistGeneObj@elementMetadata@listData$distance
}

# Calculate distances from full-length L1 to nearest gene
L1DistGene <- Dist2ClosestGene(L1CatalogGR)
L1NonZeroDistGene <- Dist2ClosestGene(L1CatalogGR_nonZero)
L1ZeroDistGene    <- Dist2ClosestGene(L1GRZero)

idxNearest <- nearest(L1CatalogGR, GRgenes, ignore.strand = T)
UpStream <- start(GRgenes)[idxNearest] < start(L1CatalogGR) 
NearestStrand <- as.vector(strand(GRgenes)[idxNearest])
table(NearestStrand, as.vector(strand(L1CatalogGR)))
chisq.test(NearestStrand, as.vector(strand(L1CatalogGR)))
table(UpStream, as.vector(strand(L1CatalogGR)))
chisq.test(NearestStrand, as.vector(strand(L1CatalogGR)))


sum(L1DistGene == 0)
hist(L1DistGene, xlab = "Distance to closest gene")
hist(L1DistGene, xlab = "Distance to closest gene", breaks = seq(0, 3*10^6, 1000))
L1DistGeneDens <- density(L1DistGene, from = 0)

# Calculate distances from full-length reference L1 to nearest gene
L1RefDistGene <- Dist2ClosestGene(L1CatalogGR_Ref)
L1RefDistGeneDens <- density(L1RefDistGene, from = 0)

# Calculate distances from fragment L1 to nearest gene
L1DistGene_Fragm <- Dist2ClosestGene(L1FragmGR)
hist(L1DistGene_Fragm, xlab = "Distance to closest gene")
hist(L1DistGene_Fragm, xlab = "Distance to closest gene", 
     breaks = seq(0, 5*10^6, 1000))
L1DistGeneDens_Fragm <- density(L1DistGene_Fragm, from = 0)

# Get the distance closest gene for full-length L1 that are not in the catalog 
L1DistGene_FullnotCat <- Dist2ClosestGene(L1RefGRFullnotCat)

# Perform regression to determine whether distance to closest gene depends on
# fragment length
FitDistVsLength <- glm(L1DistGene_Fragm ~ width(L1FragmGR))
summary(FitDistVsLength)
QRFitDistVsLength <- rq(L1DistGene_Fragm ~ width(L1FragmGR), tau = 0.75)
summary(QRFitDistVsLength)

# Calculate distances from close fragment L1 to nearest gene
L1DistGene_Fragm_close   <- Dist2ClosestGene(L1FragmGR_close)
L1DistGene_Fragm_diverge <- Dist2ClosestGene(L1FragmGR_diverge)
mean(L1DistGene_Fragm_close)
mean(L1DistGene_Fragm_diverge)

# Plot smoothed densities for catalog and fragment distance distribution
par(mfrow = c(1, 1))
plot(L1DistGeneDens$x,L1DistGeneDens$y, type = "l", ylab = "Density", 
     xlab = "Distance to closest gene", col = "blue")
lines(L1DistGeneDens_Fragm$x, L1DistGeneDens_Fragm$y, col = "red")
lines(L1RefDistGeneDens$x, L1RefDistGeneDens$y, col = "green")
legend("topright", legend = c("catalog", "reference", "fragment"), col = c("blue", "green", "red"),
       lty = c(1,1, 1))
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistDensities.pdf')

# Plot quantiles against each other
par(mfrow = c(1, 1))
qqplot(L1DistGene_Fragm, L1DistGene, ylab = "Distance full-length L1 to gene",
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

##########
#  Explore gene ontologies
##########

# SubsetGO <- function(GOs, Evidence, Ontology) {
#   IDs <- sapply(GOs, function(x) x$GOID)
#   Evs <- sapply(GOs, function(x) x$Evidence)
#   Os <- sapply(GOs, function(x) x$Ontology)
#   return(IDs[Evs %in% Evidence & Os %in% Ontology])
# }
# 
# # Get IDs of genes overlapping with catalog L1
# GeneIDs_Cat <- GRgenes@elementMetadata@listData$gene_id[blnGeneOverlapCat]
# writeLines(GeneIDs_Cat, "D:/L1polymORF/Data/GenesWithL1Insertions")
# GeneIDs_Cat_GOs <- select(org.Hs.eg.db, GeneIDs_Cat, "GO" )
# blnSubset <- GeneIDs_Cat_GOs$ONTOLOGY == "BP" & 
#   GeneIDs_Cat_GOs$EVIDENCE %in% c("IDA", "TAS")
# GeneIDs_Cat_GOs <- GeneIDs_Cat_GOs[blnSubset,]
# GOTermsBimap <-GOTERM[GeneIDs_Cat_GOs$GO]
# GoTermTable <- toTable(GOTermsBimap)
# 
# unique(GoTermTable$Term)

############################
#                          #
#     Sample random L1     #
#                          #
############################

# Determine sample size
NrSamples <- 1000

# Initialize vectors and matrices for sampled quantities
SampledGeneIntersectCount <- rep(NA, NrSamples)
SampledGeneDist           <- matrix(nrow = length(L1CatalogGR), 
                                    ncol = NrSamples)
# Create sampled ranges
for (j in 1:NrSamples) {
  SampledRanges <- lapply(which(L1CountMatched > 0), function(i) {
    Starts <- sample(ChromLengths[i], L1CountMatched[i], replace = T)
    Strands <- sample(c("+", "-"), L1CountMatched[i], replace = T)
    GRanges(seqnames = ChromNames[i], 
            ranges = IRanges(start = Starts, end = Starts + 6000),
            strand = Strands)
  })
  SampledRanges <- unlist(GRangesList(SampledRanges))
  SampledOverlapGene <- countOverlaps(SampledRanges, GRgenes) 
  SampledGeneIntersectCount[j] <- sum(SampledOverlapGene)
  SampledDistGeneObj <- distanceToNearest(SampledRanges, GRgenes) 
  SampledGeneDist[,j] <- SampledDistGeneObj@elementMetadata@listData$distance
  
}

# Plot histogram with sampled number of L1s intersecting with genes
hist(SampledGeneIntersectCount, xlab = "Number of random insertions in genes",
     main = "")
segments(sum(blnL1OverlapGene), 0, sum(blnL1OverlapGene), 500, col = "red")
sum(SampledGeneIntersectCount <= sum(blnL1OverlapGene)) / NrSamples
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntronRandom.pdf')

# Plot histogram with sampled mean distance to gene
SampledMeanDist <- colMeans(SampledGeneDist)
hist(SampledMeanDist, xlab = "Mean distance of random insertions to closest gene",
     main = "")
segments(mean(L1DistGene), 0, mean(L1DistGene), 500, col = "red")
sum(SampledMeanDist <= mean(L1DistGene)) / NrSamples
CreateDisplayPdf('D:/L1polymORF/Figures/L1MeanGeneDistRandom.pdf')

# Plot histogram with sampled median distance to gene
SampledMedDist <- apply(SampledGeneDist, 2, median)
hist(SampledMedDist, xlab = "Median distance of random insertions to closest gene",
     main = "")
segments(median(L1DistGene), 0, median(L1DistGene), 500, col = "red")
sum(SampledMedDist <= median(L1DistGene)) / NrSamples

###########################################
#                                         #
#     Sample random L1 from fragments     #
#                                         #
###########################################

# Determine sample size
NrSamples <- 1000

L1GR_actual <- L1CatalogGR
L1DistGene <- Dist2ClosestGene(L1GR_actual)

# Initialize vectors and matrices for sampled quantities
SampledGeneIntersectCount <- rep(NA, NrSamples)
SampledGeneDist           <- matrix(nrow = length(L1GR_actual), 
                                    ncol = NrSamples)
L1ToSampleFrom <- resize(L1FragmGR, width = 6064, fix = "center")

# Create sampled ranges
cat("Sampling L1 ranges from fragments\n")
for (j in 1:NrSamples) {
  SampledIndices <- sample(length(L1FragmGR), length(L1CatalogGR))
  SampledRanges  <- L1ToSampleFrom[SampledIndices]
  SampledOverlapGene <- countOverlaps(SampledRanges, GRgenes) 
  SampledGeneIntersectCount[j] <- sum(SampledOverlapGene)
  SampledGeneDist[,j] <- Dist2ClosestGene(SampledRanges)
}

# Plot histogram with sampled number of L1s intersecting with genes
hist(SampledGeneIntersectCount, xlab = "Number of sampled L1 in genes",
     breaks = seq(5, 45, 2))
segments(sum(blnL1OverlapGene), 0, sum(blnL1OverlapGene), 500, col = "red")
sum(SampledGeneIntersectCount <= sum(blnL1OverlapGene)) / NrSamples

# Plot histogram with sampled mean distance to gene
SampledMeanDist <- colMeans(SampledGeneDist)
hist(SampledMeanDist, xlab = "Mean distance of sampled L1 fragments to closest gene",
     main = "")
segments(mean(L1DistGene), 0, mean(L1DistGene), 500, col = "red")
CreateDisplayPdf('D:/L1polymORF/Figures/L1propIntronRandom.pdf')

sum(SampledMeanDist <= mean(L1DistGene)) / NrSamples
hist(apply(SampledGeneDist, 2, FUN = function(x) sum(x == 0)))

# Plot histogram with sampled maximum distance to gene
par(mfrow = c(1,1))
SampledMaxDist <- apply(SampledGeneDist, 2, max)
hist(SampledMaxDist, xlab = "Max distance of sampled L1 to genes")
segments(max(L1DistGene), 0, max(L1DistGene), 500, col = "red")
sum(SampledMaxDist <= max(L1DistGene)) / NrSamples


# Create Quantile distributions
QVect <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
QuantileMat_Sampled <- apply(SampledGeneDist, 2, 
                             FUN = function(x) quantile(x, QVect))
dim(QuantileMat_Sampled)
L1DistGeneQuantiles <- quantile(L1DistGene, QVect)
par(mfrow = c(3,2), oma = c(3,3,0,0))
for (i in 2:length(QVect)){
  hist(QuantileMat_Sampled[i, ], xlab = "", ylab = "", 
       main = paste("Quantile =", QVect[i]))
  segments(L1DistGeneQuantiles[i], 0, L1DistGeneQuantiles[i], 
           500, col = "red")
}
mtext("Distance to closest gene", side = 1, outer = T)
mtext("Frequency", side = 2, outer = T)
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQuantiles.pdf')

# Plot histogram with sampled median distance to gene
SampledMedDist <- apply(SampledGeneDist, 2, median)
hist(SampledMedDist, xlab = "Median distance of sampled L1 to genes")
segments(median(L1DistGene), 0, median(L1DistGene), 500, col = "red")
sum(SampledMedDist <= median(L1DistGene)) / NrSamples

###########################################
#                                         #
#     Sample random L1 from divergent fragments     #
#                                         #
###########################################

# Determine sample size
NrSamples <- 1000
#L1DistGene_Fragm_close   <- Dist2ClosestGene(L1FragmGR_close)
#L1DistGene_Fragm_diverge <- Dist2ClosestGene(L1FragmGR_diverge)

# Initialize vectors and matrices for sampled quantities
SampledGeneIntersectCount <- rep(NA, NrSamples)
SampledGeneDist           <- matrix(nrow = length(L1FragmGR_close), 
                                    ncol = NrSamples)
L1ToSampleFrom <- L1FragmGR_diverge

# Create sampled ranges
cat("Sampling L1 ranges from fragments\n")
for (j in 1:NrSamples) {
  SampledIndices <- sample(length(L1ToSampleFrom), length(L1FragmGR_close))
  SampledRanges  <- L1ToSampleFrom[SampledIndices]
  SampledOverlapGene <- countOverlaps(SampledRanges, GRgenes) 
  SampledGeneIntersectCount[j] <- sum(SampledOverlapGene)
  SampledGeneDist[,j] <- Dist2ClosestGene(SampledRanges)
}

# Create Quantile distributions
QVect <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
QuantileMat_Sampled <- apply(SampledGeneDist, 2, 
                             FUN = function(x) quantile(x, QVect))
L1DistGeneQuantiles <- quantile(L1DistGene_Fragm_close, QVect)
par(mfrow = c(3,2), oma = c(3,3,0,0))
for (i in 2:length(QVect)){
  hist(QuantileMat_Sampled[i, ], xlab = "", ylab = "", 
       main = paste("Quantile =", QVect[i]))
  segments(L1DistGeneQuantiles[i], 0, L1DistGeneQuantiles[i], 
           500, col = "red")
}
mtext("Distance to closest gene", side = 1, outer = T)
mtext("Frequency", side = 2, outer = T)
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistQuantilesCloseFragm.pdf')

# Plot histogram with sampled median distance to gene
SampledMedDist <- apply(SampledGeneDist, 2, median)
hist(SampledMedDist, xlab = "Median distance of sampled L1 to genes")
segments(median(L1DistGene), 0, median(L1DistGene), 500, col = "red")
sum(SampledMedDist <= median(L1DistGene)) / NrSamples



