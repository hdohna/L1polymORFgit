# The following script reads in a Hi-C data and determines how much
# different L1 insertions interact with other genomic regions

library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(irlba)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Source start script
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
Resolution <- '10kb'
WindowL <- as.numeric(strsplit(Resolution, 'kb')[[1]][1]) * 1000
NRows   <- 10^7

# Specify patterns 
NormPattern <- "VCnorm"
ExpPattern  <- "VCexpected"
NormPattern <- "KRnorm"
ExpPattern  <- "KRexpected"
MAPQFolder  <- "MAPQGE30"

#  Set paths
RepeatTablePath <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/repeatsHg19_L1HS.csv"
ChainFile38To19 <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg38ToHg19.over.chain"
ChainFile19To18 <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/hg19ToHg18.over.chain"
L1CatalogPath   <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"
HiCFolderPath   <- paste("/srv/gsfs0/projects/levinson/hzudohna/HiCData/GM12878_combined/", 
                         Resolution, "_resolution_intrachromosomal/", sep = "")
ChromLenghtPath <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/ChromLengthsHg19.Rdata"
OutputPath <- paste("/srv/gsfs0/projects/levinson/hzudohna/HiCData/GM12878_combined/", Resolution, 
                    NormPattern, "L1summary.Rdata", sep = "_")


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

# Initialize data frame that collects HiC Data
HicByL1Type <- data.frame()

# Set chromosome and window length
for (Chrom in names(ChromLengthsHg19)[names(ChromLengthsHg19) != "chrY"]) {
  cat("\n**********    Processing chromosome", Chrom, "    *************\n")
  CurrentFolder <- paste(HiCFolderPath, Chrom, MAPQFolder, sep = "/")
  
  # Create a genomic ranges of all HiC windows
  WStarts <- seq(0, ChromLengthsHg19[Chrom], WindowL)
  HiCGR   <- GRanges(seqnames = Chrom, 
                     ranges = IRanges(start = WStarts, width = WindowL))
  
  overlapL1Fragm <- findOverlaps(L1GRhg19_fragm, HiCGR)
  overlapL1Cat   <- findOverlaps(L1GRhg19_cat, HiCGR)
  idxHicGR       <- c(overlapL1Fragm@to, overlapL1Cat@to)
  WStartsL1      <- WStarts[idxHicGR]
  cat(length(WStartsL1), "ranges intersect with L1s\n")
  
  # Get lists of files 
  RawMatFile <- list.files(CurrentFolder, full.names = T, 
                           pattern = paste(Resolution, "RAWobserved", sep = "."))
  StVFile   <- list.files(CurrentFolder, full.names = T, 
                          pattern = paste(Resolution, NormPattern, sep = "."))
  ExpFile   <- list.files(CurrentFolder, full.names = T, 
                          pattern = paste(Resolution, ExpPattern, sep = "."))
  
  # Read in Hi-C data and standardization vectors
  StVect  <- read.table(StVFile)
  ExpVect <- read.table(ExpFile)
  
  # Create a standardization vector for the number of genes per range
  PromGR    <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
  PromCount <- countOverlaps(HiCGR, PromGR)
  PromCount <- PromCount[-1]
  cat("Number of HiC genomic ranges:", length(HiCGR), "\n")
  cat("Length of norm vector:", nrow(StVect), "\n")
  
  HiCAgg       <- data.frame()
  Lines2Skip   <- 0 
  LinesRead    <- NRows
  TotLinesRead <- 0
  MinLeft1     <- Inf
  MaxLeft1     <- 0
  MinLeft2     <- Inf
  MaxLeft2     <- 0
  while(LinesRead == NRows){
    
    # Read in matrix of H-C values and standardization vectors
    HiCMat           <- read.table(RawMatFile, nrows = NRows, skip = Lines2Skip)
    LinesRead  <- nrow(HiCMat)
    colnames(HiCMat) <- c("Left1", "Left2", "RawReads")
    blnLeft1InL1 <- HiCMat$Left1 %in% WStartsL1 
    blnLeft2InL1 <- HiCMat$Left2 %in% WStartsL1 
    MinLeft1 <- min(MinLeft1, HiCMat$Left1)
    MaxLeft1 <- max(MaxLeft1, HiCMat$Left1)
    MinLeft2 <- min(MinLeft2, HiCMat$Left2)
    MaxLeft2 <- max(MaxLeft2, HiCMat$Left2)
    
    HiCAggNew1 <- data.frame()
    HiCAggNew2 <- data.frame()
    
    if (any(blnLeft1InL1)) {
      
      # Subset and standardize Hi-C data
      HiCMat1 <- HiCMat[blnLeft1InL1, ]
      idx1             <- HiCMat1[,1] / WindowL + 1
      idx2             <- HiCMat1[,2] / WindowL + 1
      idx3             <- (HiCMat1[,2] - HiCMat1[,1]) / WindowL + 1
      HiCMat1$NormReads <- HiCMat1[,3] / StVect[idx1,] / StVect[idx2,] 
      HiCMat1$NormEO    <- HiCMat1$NormReads / ExpVect[idx3,] 
      HiCMat1$NormProm1 <- HiCMat1$NormEO * PromCount[idx1]
      HiCMat1$NormProm2 <- HiCMat1$NormEO * PromCount[idx2]
      HiCAggNew1 <- aggregate(cbind(NormEO, NormProm2) ~ Left1, data = HiCMat1, FUN = sum)
    }
    if (any(blnLeft2InL1)) {
      
      # Subset and standardize Hi-C data
      HiCMat2 <- HiCMat[blnLeft2InL1, ]
      idx1             <- HiCMat2[,1] / WindowL + 1
      idx2             <- HiCMat2[,2] / WindowL + 1
      idx3             <- (HiCMat2[,2] - HiCMat2[,1]) / WindowL + 1
      HiCMat2$NormReads <- HiCMat2[,3] / StVect[idx1,] / StVect[idx2,] 
      HiCMat2$NormEO    <- HiCMat2$NormReads / ExpVect[idx3,] 
      HiCMat2$NormProm1 <- HiCMat2$NormEO * PromCount[idx1]
      HiCMat2$NormProm2 <- HiCMat2$NormEO * PromCount[idx2]
      HiCAggNew2 <- aggregate(cbind(NormEO, NormProm1) ~ Left2, data = HiCMat2, FUN = sum)
      colnames(HiCAggNew2)[colnames(HiCAggNew2) == "Left2"] <- "Left1"
      colnames(HiCAggNew2)[colnames(HiCAggNew2) == "NormProm1"] <- "NormProm2"
    }   
      # Aggregate interchromosomal interaction
      HiCAgg    <- rbind(HiCAgg, HiCAggNew1, HiCAggNew2)
      HiCAgg    <- aggregate(cbind(NormEO, NormProm2) ~ Left1, data = HiCAgg, FUN = sum)

    # Update iteration variables and produce status message
    cat("Processed line", Lines2Skip + 1, "to", Lines2Skip + nrow(HiCMat), "\n")
    Lines2Skip <- Lines2Skip + LinesRead 
    TotLinesRead <- TotLinesRead + LinesRead
    cat("Total lines read is", TotLinesRead, "\n")
    cat("Current HiCAgg rows:", nrow(HiCAgg), "\n")
  }
  
  # Report minimum and maximum Left1
  cat("Minimum Left1:", MinLeft1, "\n")
  cat("Maximum Left1:", MaxLeft1, "\n")
  cat("Minimum Left2:", MinLeft2, "\n")
  cat("Maximum Left2:", MaxLeft2, "\n")
  cat("Chromosome length:", MaxLeft2, "\n")
  
  
  # Aggregate total HiC interaction value by L1
  StartMatch <- match(HiCAgg$Left1, start(HiCGR))
  HicByL1TypeNew1 <- AggregateValsBy2GRangesSet(L1GRhg19_cat, L1GRhg19_fragm, 
                      HiCGR[StartMatch], HiCAgg$NormProm2, 
                      Type12Names = c("full", "fragm"), ValueName = "NormPromSum", 
                      TypeName = "L1type")
  HicByL1TypeNew2 <- AggregateValsBy2GRangesSet(L1GRhg19_cat, L1GRhg19_fragm, 
                      HiCGR[StartMatch], HiCAgg$NormEO, 
                      Type12Names = c("full", "fragm"), ValueName = "NormEOsum", 
                      TypeName = "L1type")
  HicByL1TypeNew <- merge(HicByL1TypeNew1, HicByL1TypeNew2)
  HicByL1Type    <- rbind(HicByL1Type, HicByL1TypeNew)
  cat("Current HicByL1Type rows:", nrow(HicByL1Type), "\n")
}

############
#  Test for difference in aggregation level 
############

# Overall test
t.test(NormEOsum   ~ L1type, HicByL1Type)
t.test(NormPromSum ~ L1type, HicByL1Type)

# Test per chromosome
pValPerChrom <- sapply(names(ChromLengthsHg19), function(x){
  HicSubset <- HicByL1Type[HicByL1Type$chromosome == x, ]
  L1Count <- table(HicSubset$L1type)
  if (min(L1Count) > 3 & length(L1Count) == 2){
    c(pEO = t.test(NormEOsum   ~ L1type, HicSubset)$p.value,
      pProm = t.test(NormPromSum ~ L1type, HicSubset)$p.value)
  } else {
    c(pEO = NA, pProm = NA)
  }
})

HiCSumPerChrom <- aggregate(cbind(NormPromSum, NormEOsum) ~ L1type + chromosome, HicByL1Type, FUN = mean)

# Save workspace
save.image(OutputPath)
