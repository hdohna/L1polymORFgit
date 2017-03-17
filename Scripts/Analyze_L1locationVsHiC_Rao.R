# The following script reads in a Hi-C data and determines how much
# different L1 insertions interact with other genomic regions

library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(irlba)


############
#  Set parameters
############

#  Set paths
RepeatTablePath <- "D:/L1polymORF/Data/repeatsHg38_L1HS.csv"
ChainFile38To19 <- "D:/L1polymORF/Data/hg38ToHg19.over.chain"
ChainFile19To18 <- "D:/L1polymORF/Data/hg19ToHg18.over.chain"
L1CatalogPath   <- "D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"
HiCFolderPath   <- "D:/L1polymORF/Data/HiCData"

# Set chromosome and window length
Chrom   <- "chr1"
WindowL <- 50000

############
#  Process L1 data
############

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv(RepeatTablePath)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
L1RefGR_hg19 <- liftOver(L1RefGR, 
                         chain = import.chain(ChainFile38To19))
NrRanges <- sapply(L1RefGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1RefGR_hg19 <- unlist(L1RefGR_hg19[idxUniqueMapped])
L1GRhg19_fragm <- L1RefGR_hg19[width(L1RefGR_hg19) < 5900]


# Path to L1 catalogue file 
L1Catalogue <- read.csv(L1CatalogPath, as.is = T)


L1Catalogue$Allele[is.na(L1Catalogue$Allele)] <- 1

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef1         <- (L1Catalogue$end_HG38 - L1Catalogue$start_HG38) > 6000 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Create genomic ranges for catalog L1
L1CatalogGR <- GRanges(seqnames = L1CatalogL1Mapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogL1Mapped$start_HG38,
                                                     L1CatalogL1Mapped$end_HG38),
                                        end = pmax(L1CatalogL1Mapped$start_HG38,
                                                   L1CatalogL1Mapped$end_HG38)),
                       strand = L1CatalogL1Mapped$strand_L1toRef)
L1CatGR_hg19 <- liftOver(L1CatalogGR, 
                         chain = import.chain(ChainFile38To19))
NrRanges <- sapply(L1CatGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1CatGR_hg19 <- unlist(L1CatGR_hg19[idxUniqueMapped])


############
#  Process Hi-C data
############

# Get lists of files 
RawMatFile <- list.files(HiCFolderPath, full.names = T, 
                         pattern = "RAWobserved")
StVFile   <- list.files(HiCFolderPath, full.names = T, pattern = "kb.VCnorm")
ExpFile   <- list.files(HiCFolderPath, full.names = T, pattern = ".VCexpected")

# Read in Hi-C data and standardization vectors
HiCMat           <- read.table(RawMatFile)
colnames(HiCMat) <- c("Left1", "Left2", "RawReads")
StVect           <- read.table(StVFile)
ExpVect          <- read.table(ExpFile)
dim(HiCMat)

# Standardize Hi-C data and standardization vectors
idx1             <- HiCMat[,1] / WindowL + 1
idx2             <- HiCMat[,2] / WindowL + 1
idx3             <- (HiCMat[,2] - HiCMat[,1]) / WindowL + 1
HiCMat$NormReads <- HiCMat[,3] / StVect[idx1,] / StVect[idx2,]
HiCMat$NormEO    <- HiCMat$NormReads / ExpVect[idx3,]

# Aggregate interchromosomal interaction
HiCAgg <- aggregate(NormEO ~ Left1, data = HiCMat, FUN = sum)
min(HiCMat$Left1)
HiCGR <- GRanges(seqnames = Chrom, 
                 IRanges(start = HiCAgg$Left1,
                 width = WindowL))

############
#  Intersect HiC data with L1
############

# Aggregate total HiC interaction value by L1
HicByL1Type <- AggregateValsBy2GRangesSet(L1CatGR_hg19, L1GRhg19_fragm, HiCGR, 
   HiCAgg$NormEO, Type12Names = c("full", "fragm"), ValueName = "NormEOsum", 
   TypeName = "L1type")
t.test(NormEOsum ~ L1type, HicByL1Type)
boxplot(NormEOsum ~ L1type, HicByL1Type)

############
#  Perform PCA
############

# # Create a sparse matrix of contact
# M <- sparseMatrix(i = idx1, j = idx2, x = HiCMat$NormEO, symmetric = T)
# SVDM <- irlba(M, nv = 1, nu = 0, center = colMeans(M), fastpath = FALSE)
# PC1 <- M %*% SVDM$v
# length(PC1)
# plot(PC1, ylim = c(-0.1, 3))
# plot(PC1[1:200], ylim = c(-0.1, 3), type = "l")
# 
# # Turn Hi-C data into a GRanges object
# idxV  <- 0:nrow(StVect)* WindowL
# HiCGR <- GRanges(seqnames = Chrom, IRanges(start = idxV[-length(idxV)],
#                                            end = idxV[-1]))
# # Find overlap with fragments and catalog elements
# L1RefGR_hg19_fragm <- L1RefGR_hg19[width(L1RefGR_hg19) < 5000]
# idxOverlap_fragm   <- findOverlaps(HiCGR, L1RefGR_hg19_fragm)
# idxOverlap_Cat     <- findOverlaps(HiCGR, L1CatGR_hg19)

