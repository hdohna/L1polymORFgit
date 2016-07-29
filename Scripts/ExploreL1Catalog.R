# The script below explore a catalog of full-length L1

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(ShortRead)
# library(csaw)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"

# Length of Flanking sequence to be used for alignment
FlankLength <- 200

# Maximum fragment length
MaxFragLength <- 3000


######################################
#                                    #
#    Read & process L1 catalog       #
#                                    #
######################################

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS",]
RepeatTable <- RepeatTable[abs(RepeatTable$genoEnd - RepeatTable$genoStart) <=
                             MaxFragLength,]

# Create genomic ranges for L1 fragments
L1FragmGR <- GRanges(seqnames = RepeatTable$genoName,
                     ranges = IRanges(start = RepeatTable$genoStart,
                                      end = RepeatTable$genoEnd),
                     strand = RepeatTable$strand)
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
L1OverlapGene <- countOverlaps(L1CatalogGR, 
                               genes(TxDb.Hsapiens.UCSC.hg38.knownGene)) 
L1OverlapPromoter <- countOverlaps(L1CatalogGR, 
   promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 5000)) 
sum(L1OverlapExon)
sum(L1OverlapGene)
sum(L1OverlapPromoter)

##########
#  Calculate distances to genes
##########

# Count overlaps with exons and promoters
L1DistGeneObj <- distanceToNearest(L1CatalogGR, 
                                genes(TxDb.Hsapiens.UCSC.hg38.knownGene)) 
L1DistGene <- L1DistGeneObj@elementMetadata@listData$distance
sum(L1DistGene == 0)
hist(L1DistGene, xlab = "Distance to closest gene")

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

