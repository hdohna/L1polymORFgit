# The following script reads in a Hi-C matrix and determines well-connected stretchs

library(GenomicRanges)
library(rtracklayer)

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")

# Create genomic ranges for L1 fragments, match them to distances to get distance
# to consensus per fragment
L1RefGR <- GRanges(seqnames = RepeatTable$genoName,
                   ranges = IRanges(start = RepeatTable$genoStart,
                                    end = RepeatTable$genoEnd),
                   strand = RepeatTable$strand)
ChainFile38To19   <- "D:/L1polymORF/Data/hg38ToHg19.over.chain"
ChainFile19To18   <- "D:/L1polymORF/Data/hg19ToHg18.over.chain"
L1RefGR_hg19 <- liftOver(L1RefGR, 
                         chain = import.chain(ChainFile38To19))
NrRanges <- sapply(L1RefGR_hg19, length)
idxUniqueMapped <- NrRanges == 1
L1RefGR_hg19 <- unlist(L1RefGR_hg19[idxUniqueMapped])
L1RefGR_hg18 <- liftOver(L1RefGR_hg19, 
                         chain = import.chain(ChainFile19To18))
NrRanges <- sapply(L1RefGR_hg18, length)
idxUniqueMapped <- NrRanges == 1
L1RefGR_hg18 <- unlist(L1RefGR_hg18[idxUniqueMapped])


# Path to L1 catalogue file 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)


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
L1CatGR_hg18 <- liftOver(L1CatGR_hg19, 
                         chain = import.chain(ChainFile19To18))
NrRanges <- sapply(L1CatGR_hg18, length)
idxUniqueMapped <- NrRanges == 1
L1CatGR_hg18 <- unlist(L1CatGR_hg18[idxUniqueMapped])


# Get lists of files 
EijFiles <- list.files("D:/L1polymORF/Data/hES/nij/", full.names = T)
#EijFiles <- list.files("D:/L1polymORF/Data/hES/uij/", full.names = T)
HiCData <- data.frame()
#EijFile <- EijFiles[1]

# Specify the number of neighbors to consider
NrNeighbors <- 4
for (EijFile in EijFiles){
  
  # Read in Eij
  cat("Processing", EijFile, "\n")
  Eij <- read.table(EijFile)
  
  # Mask
  NrBins  <- nrow(Eij)
  idxVect <- 1:NrBins
  MasKMat <- sapply(idxVect, function(i) {
    idxVect %in% (max(1, i - NrNeighbors):min(NrBins, i + NrNeighbors))
  })
  
  Eij <- rowSums(Eij * MasKMat)
  
  # Create new data 
  NewData <- data.frame(
    chrom = strsplit(EijFile, "\\.")[[1]][2],
    start = seq(1, length(Eij) * 40000, 40000),
    end   = seq(40000, length(Eij) * 40000, 40000),
    Eij   = Eij)

  # Append the new data to the Hi-C table
  HiCData <- rbind(HiCData, NewData)
}

# Turn Hi-C data into a GRanges object
HiCGR <- makeGRangesFromDataFrame(HiCData)

# Find overlap with fragments and catalog elements
L1RefGR_hg18_fragm <- L1RefGR_hg18[width(L1RefGR_hg18) < 5000]
idxOverlap_fragm   <- findOverlaps(HiCGR, L1RefGR_hg18_fragm)
idxOverlap_Cat     <- findOverlaps(HiCGR, L1CatGR_hg18)

# Group HiC data by overlap
HiCByOverlap <- data.frame(HiC = c(HiCData$Eij[idxOverlap_fragm@queryHits],
                                   HiCData$Eij[idxOverlap_Cat@queryHits]),
                           L1Type = c(rep("Fragm", length(idxOverlap_fragm@queryHits)),
                                      rep("Full", length(idxOverlap_Cat@queryHits))))
boxplot(HiCByOverlap$HiC ~ HiCByOverlap$L1Type)

# Distinguish full-length and fragment
idxOverlap <- findOverlaps(HiCGR, L1RefGR_hg18)
idxOverlap@subjectHits

blnFullL1 <- width(L1RefGR_hg18[idxOverlap@subjectHits]) > 6000
boxplot(HiCData$Eij[idxOverlap@queryHits] ~ blnFullL1)
aggPerWindow <- aggregate(blnFullL1 ~ idxOverlap@queryHits, FUN = mean)
aggPerWindow$Eij <- HiCData$Eij[aggPerWindow$`idxOverlap@queryHits`]
plot(aggPerWindow$Eij, aggPerWindow$blnFullL1)
