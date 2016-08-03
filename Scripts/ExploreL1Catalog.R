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

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"

# Path to L1 alignment file 
AlignFileName <- 'D:/L1polymORF/Data/L1HSSequences_L100_withConsens_aligned.fas'

# Length of Flanking sequence to be used for alignment
FlankLength <- 200

# Maximum fragment length
MaxFragLength <- 5000


######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

# Read in alignment
L1HSAlign    <- read.fasta(AlignFileName)

# Specify motif ttagtgggtg (nucleotide sequence after ACA) and find ACA 
# locus in sequences
MotifAfterACA <- c("t", "t", "a", "g", "t", "g", "g", "g", "t", "g")
MotifL <- length(MotifAfterACA)
L1ACALocus   <- sapply(L1HSAlign, function(x) {
  SeqCollapsed <- x[x != "-"]
  idxACA <- which(sapply(1:length(SeqCollapsed), function(i){
    all(SeqCollapsed[i:(i + MotifL - 1)] == MotifAfterACA)
  })) - 3
  if (length(idxACA) == 1){
    paste(SeqCollapsed[idxACA:(idxACA + 2)], collapse = "")
  } else {
    NA
  }
})
L1HSAlign    <- read.dna(AlignFileName, format = "fasta")
Dist2Consens  <- dist.dna(L1HSAlign, as.matrix = T, pairwise.deletion = T)
Dist2Consens  <- Dist2Consens[1,]
# L1HSAlign    <- read.fasta(AlignFileName)
# FragmLength <- sapply(L1HSAlign, function(x) sum(x != "-"))

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS",]
RepeatTable <- RepeatTable[abs(RepeatTable$genoEnd - RepeatTable$genoStart) <=
                             MaxFragLength,]

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1FragmGR <- GRanges(seqnames = RepeatTable$genoName,
                     ranges = IRanges(start = RepeatTable$genoStart,
                                      end = RepeatTable$genoEnd),
                     strand = RepeatTable$strand)
L1FragmNames <- paste(seqnames(L1FragmGR), start(L1FragmGR), end(L1FragmGR),
                      strand(L1FragmGR), sep = "_")
NameMatch <- match(L1FragmNames, names(Dist2Consens))
Dist2Consens <- Dist2Consens[NameMatch]
hist(Dist2Consens, breaks = seq(0, 2, 0.05))
plot(width(L1FragmGR), Dist2Consens)
sum(Dist2Consens < 0.4, na.rm = T)

# Get 2 sets of fragments: close to the L1HS consensus and divergent
L1FragmGR_close   <- L1FragmGR[which(Dist2Consens < 0.5)]
L1FragmGR_diverge <- L1FragmGR[which(Dist2Consens >= 0.5)]

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
   ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                 L1CatalogL1Mapped$end_HG38),
                    end = pmax(L1CatalogL1Mapped$start_HG38,
                               L1CatalogL1Mapped$end_HG38)),
                      strand = L1CatalogL1Mapped$Strand)


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
ChromNames <- paste("chr", c(1:22, "X", "Y"), sep = "")
ChromLengths <- sapply(ChromNames, function(x) length(BSgenome.Hsapiens.UCSC.hg38[[x]]))

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
plot(ChromLengths, L1CountMatched, xlab = "Chromosome length", ylab = "L1 count")
abline(coefficients(PoissonFit))

# Plot chromosome length vs. L1 fragment count
plot(ChromLengths, L1CountMatched_Fragm, xlab = "Chromosome length", 
     ylab = "L1 fragment count")

##########
#  Count L1 intersecting with exons, genes and propmoters
##########

# Get ranges of genes
GRgenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Count overlaps with exons and promoters
L1OverlapExon <- countOverlaps(L1CatalogGR, 
                               exons(TxDb.Hsapiens.UCSC.hg38.knownGene)) 
L1OverlapGene <- countOverlaps(L1CatalogGR, GRgenes)
L1OverlapPromoter <- countOverlaps(L1CatalogGR, 
   promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 5000)) 
sum(L1OverlapExon)
sum(L1OverlapGene)
sum(L1OverlapPromoter)
mean(L1OverlapGene)

# Perform logistic regression to determine whether fragment length predicts
# probability of gene overlap
L1OverlapGene_Fragm <- countOverlaps(L1FragmGR, GRgenes)
L1OverlapGene_Fragm <- L1OverlapGene_Fragm > 0
LogRegFit <- glm(L1OverlapGene_Fragm ~ width(L1FragmGR),
                 family = binomial)
summary(LogRegFit)
WidthOrder <- order(width(L1FragmGR), decreasing = T)
plot(width(L1FragmGR)[WidthOrder], L1OverlapGene_Fragm[WidthOrder])
lines(width(L1FragmGR)[WidthOrder], 
     fitted.values(LogRegFit)[WidthOrder], type = "l")

# Get proportion overlap 
WidthCut <- cut(width(L1FragmGR), breaks = seq(0, 5000, 500))
PropInGene <- aggregate(L1OverlapGene_Fragm ~ WidthCut, FUN = mean)
AvWidth    <- aggregate(width(L1FragmGR) ~ WidthCut, FUN = mean)
plot(AvWidth$`width(L1FragmGR)`, PropInGene$L1OverlapGene_Fragm,
     xlab = "Fragment size [bp]", ylab = "Proportion in exons")
lines(width(L1FragmGR)[WidthOrder], fitted.values(LogRegFit)[WidthOrder])
segments(0, mean(L1OverlapGene), 10^6, mean(L1OverlapGene), col = "red",
         lty = 2)
sum(width(L1FragmGR) > 3000)

##########
#  Calculate distances to genes
##########

# Auxiliary function to get distances to closest gene
Dist2ClosestGene <- function(GR){
  DistGeneObj <- distanceToNearest(GR, GRgenes) 
  DistGeneObj@elementMetadata@listData$distance
}

# Calculate distances from full-length L1 to nearest gene
L1DistGene <- Dist2ClosestGene(L1CatalogGR)
sum(L1DistGene == 0)
hist(L1DistGene, xlab = "Distance to closest gene")
hist(L1DistGene, xlab = "Distance to closest gene", breaks = seq(0, 3*10^6, 1000))
L1DistGeneDens <- density(L1DistGene, from = 0)

# Calculate distances from fragment L1 to nearest gene
L1DistGene_Fragm <- Dist2ClosestGene(L1FragmGR)
hist(L1DistGene_Fragm, xlab = "Distance to closest gene")
hist(L1DistGene_Fragm, xlab = "Distance to closest gene", breaks = seq(0, 5*10^6, 1000))
L1DistGeneDens_Fragm <- density(L1DistGene_Fragm, from = 0)

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
legend("topright", legend = c("catalog", "fragment"), col = c("blue", "red"),
       lty = c(1,1))
CreateDisplayPdf('D:/L1polymORF/Figures/L1geneDistDensities.pdf')

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
hist(SampledGeneIntersectCount, xlab = "Number of sampled L1 in genes")
segments(sum(L1OverlapGene), 0, sum(L1OverlapGene), 500, col = "red")
sum(SampledGeneIntersectCount <= sum(L1OverlapGene)) / NrSamples

# Plot histogram with sampled mean distance to gene
SampledMeanDist <- colMeans(SampledGeneDist)
hist(SampledMeanDist, xlab = "Mean distance of sampled L1 to genes")
segments(mean(L1DistGene), 0, mean(L1DistGene), 500, col = "red")
sum(SampledMeanDist <= mean(L1DistGene)) / NrSamples

# Plot histogram with sampled median distance to gene
SampledMedDist <- apply(SampledGeneDist, 2, median)
hist(SampledMedDist, xlab = "Median distance of sampled L1 to genes")
segments(median(L1DistGene), 0, median(L1DistGene), 500, col = "red")
sum(SampledMedDist <= median(L1DistGene)) / NrSamples


###########################################
#                                         #
#     Sample random L1 from fragments     #
#                                         #
###########################################

# Determine sample size
NrSamples <- 1000

# Initialize vectors and matrices for sampled quantities
SampledGeneIntersectCount <- rep(NA, NrSamples)
SampledGeneDist           <- matrix(nrow = length(L1CatalogGR), 
                                    ncol = NrSamples)
L1ToSampleFrom <- resize(L1FragmGR, width = 6064, fix = "center")

# Create sampled ranges
cat("Sampling L1 ranges from fragments\n")
for (j in 1:NrSamples) {
  SampledIndices <- sample(length(L1FragmGR), length(L1CatalogGR))
  SampledRanges  <- L1ToSampleFrom[SampledIndices]
  SampledOverlapGene <- countOverlaps(SampledRanges, GRgenes) 
  SampledGeneIntersectCount[j] <- sum(SampledOverlapGene)
  SampledDistGeneObj  <- distanceToNearest(SampledRanges, GRgenes) 
  SampledGeneDist[,j] <- SampledDistGeneObj@elementMetadata@listData$distance
}

# Plot histogram with sampled number of L1s intersecting with genes
hist(SampledGeneIntersectCount, xlab = "Number of sampled L1 in genes",
     breaks = seq(5, 45, 2))
segments(sum(L1OverlapGene), 0, sum(L1OverlapGene), 500, col = "red")
sum(SampledGeneIntersectCount <= sum(L1OverlapGene)) / NrSamples

# Plot histogram with sampled mean distance to gene
SampledMeanDist <- colMeans(SampledGeneDist)
hist(SampledMeanDist, xlab = "Mean distance of sampled L1 to genes")
segments(mean(L1DistGene), 0, mean(L1DistGene), 500, col = "red")
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



