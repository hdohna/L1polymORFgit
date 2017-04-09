# The following script reads in a Hi-C data and determines how much
# different L1 insertions interact with other genomic regions

library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(irlba)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Source start script
#source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

############
#  Set parameters
############

#  Set paths
RepeatTablePath <- "D:/L1polymORF/Data/repeatsHg19_L1HS.csv"
ChainFile38To19 <- "D:/L1polymORF/Data/hg38ToHg19.over.chain"
ChainFile19To18 <- "D:/L1polymORF/Data/hg19ToHg18.over.chain"
L1CatalogPath   <- "D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"
HiCFolderPath   <- "D:/L1polymORF/Data/HiCData"
ChromLenghtPath <- "D:/L1polymORF/Data/ChromLengthsHg19.RData"

# Set window length and number of rows to be read in
WindowLchar <- '50kb'
WindowL <- as.numeric(strsplit(WindowLchar, 'kb')[[1]][1]) * 1000
NRows   <- 10^7

#  Set paths
# RepeatTablePath <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/repeatsHg19_L1HS.csv"
# ChainFile38To19 <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg38ToHg19.over.chain"
# ChainFile19To18 <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg19ToHg18.over.chain"
# L1CatalogPath   <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"
# HiCFolderPath   <- "/srv/gsfs0/projects/levinson/hzudohna/HiCData/GM12878_combined/50kb_resolution_intrachromosomal/"
# MAPQFolder      <- "MAPQGE30"
# ChromLenghtPath <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/ChromLengthsHg19.Rdata"


############
#  Process L1 data
############

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv(RepeatTablePath)

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR_hg19 <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
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
L1GRhg19_cat    <- liftOver(L1CatalogGR, 
                         chain = import.chain(ChainFile38To19))
NrRanges        <- sapply(L1GRhg19_cat, length)
idxUniqueMapped <- NrRanges == 1
L1GRhg19_cat    <- unlist(L1GRhg19_cat[idxUniqueMapped])

# Load data with chromosome length
load(ChromLenghtPath)
ChromLengthsHg19

############
#  Process Hi-C data
############

# Set chromosome and window length
Chrom   <- "chr1"
CurrentFolder <- HiCFolderPath

# Create a genomic ranges of all HiC windows
WStarts <- seq(0, ChromLengthsHg19[Chrom], WindowL)
HiCGR   <- GRanges(seqnames = Chrom, 
                   ranges = IRanges(start = WStarts, width = WindowL))

overlapL1Fragm <- findOverlaps(L1GRhg19_fragm, HiCGR)
overlapL1Cat   <- findOverlaps(L1GRhg19_cat, HiCGR)
idxHicGR       <- c(overlapL1Fragm@to, overlapL1Cat@to)
WStartsL1      <- WStarts[idxHicGR]

# Get lists of files 
RawMatFile <- list.files(CurrentFolder, full.names = T, 
                pattern = paste(WindowLchar, "RAWobserved", sep = "."))
StVFile   <- list.files(CurrentFolder, full.names = T, 
                pattern = paste(WindowLchar, "VCnorm", sep = "."))
ExpFile   <- list.files(CurrentFolder, full.names = T, 
                pattern = paste(WindowLchar, "VCexpected", sep = "."))

# Read in Hi-C data and standardization vectors
StVect  <- read.table(StVFile)
ExpVect <- read.table(ExpFile)


# Create a standardization vector for the number of genes per range
PromGR    <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
PromCount <- countOverlaps(HiCGR, PromGR)
PromCount <- PromCount[-1]

HiCAgg       <- data.frame()
Lines2Skip   <- 0 
LinesRead    <- NRows
TotLinesRead <- 0
while(LinesRead == NRows){
  
  # Read in matrix of H-C values and standardization vectors
  HiCMat           <- read.table(RawMatFile, nrows = NRows, skip = Lines2Skip)
  LinesRead  <- nrow(HiCMat)
  colnames(HiCMat) <- c("Left1", "Left2", "RawReads")
  blnInL1 <- (HiCMat$Left1 %in% WStartsL1) | (HiCMat$Left2 %in% WStartsL1) 
  HiCMat  <- HiCMat[blnInL1,]
  
  if (nrow(HiCMat) > 0) {
    # Standardize Hi-C data and standardization vectors
    idx1             <- HiCMat[,1] / WindowL + 1
    idx2             <- HiCMat[,2] / WindowL + 1
    idx3             <- (HiCMat[,2] - HiCMat[,1]) / WindowL + 1
    HiCMat$NormReads <- HiCMat[,3] / StVect[idx1,] / StVect[idx2,] 
    HiCMat$NormEO    <- HiCMat$NormReads / ExpVect[idx3,] 
    HiCMat$NormProm1 <- HiCMat$NormEO * PromCount[idx1]
    HiCMat$NormProm2 <- HiCMat$NormEO * PromCount[idx2]
    
    # Turn HiC data into a sparse matrix
    HiCSpM <- sparseMatrix(idx1, idx2)
    
    # Aggregate interchromosomal interaction
    HiCAggNew1 <- aggregate(cbind(NormEO, NormProm2) ~ Left1, data = HiCMat, FUN = sum)
    HiCAggNew2 <- aggregate(cbind(NormEO, NormProm1) ~ Left2, data = HiCMat, FUN = sum)
    colnames(HiCAggNew2)[colnames(HiCAggNew2) == "Left2"] <- "Left1"
    colnames(HiCAggNew2)[colnames(HiCAggNew2) == "NormProm1"] <- "NormProm2"
    HiCAgg    <- rbind(HiCAgg, HiCAggNew1, HiCAggNew2)
    HiCAgg    <- aggregate(cbind(NormEO, NormProm2) ~ Left1, data = HiCAgg, FUN = sum)
  }
  
  
  # Update iteration variables and produce status message
  cat("Processed line", Lines2Skip + 1, "to", Lines2Skip + nrow(HiCMat), "\n")
  Lines2Skip <- Lines2Skip + LinesRead 
  TotLinesRead <- TotLinesRead + LinesRead
  cat("Total lines read is", TotLinesRead, "\n")
}


############
#  Intersect HiC data with L1
############

# Aggregate total HiC interaction value by L1
blnHicGRAgg <- start(HiCGR) %in% HiCAgg$Left1
HicByL1Type <- AggregateValsBy2GRangesSet(L1GRhg19_cat, L1GRhg19_fragm, 
                                          HiCGR[blnHicGRAgg], 
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
# idxOverlap_Cat     <- findOverlaps(HiCGR, L1GRhg19_cat)

