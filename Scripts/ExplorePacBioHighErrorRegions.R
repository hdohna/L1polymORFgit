# This script reads in bam files and determines correlation between phred scores of 
# different reads 

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Minimum proportion of SNP in reads to be called
MinPropCall <- 0.6

# Paths (SCG4)
PathStart          <- '/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R'
PathSNPRanges      <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/SNPinRangesNA12878.RData'
PathBam_Normal_BWA <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
PathOutput         <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioHighErrorRegions.RData'

# Source start script
source(PathStart)

###############################################
#                                             #
#    Load ranges where reads were mapped      #
#                                             #
###############################################

cat("\n**********    Loading data    *************\n")

# Load data with ranges where reads were mapped and SNPs (generated in script
# PrepareSNPsNA12878.R)
load(PathSNPRanges)

# Subset GRUnion to get the ones containing SNPs 
GRUnion_withSNP <- GRUnion_withSNP[width(GRUnion_withSNP) >= 1000]

# Get reference sequence for each range
RefSeqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRUnion_withSNP)

# Get reads for different datasets
scanPars        <- ScanBamParam(which = GRUnion_withSNP, what = scanBamWhat())
ReadList_Normal_BWA <- scanBam(PathBam_Normal_BWA, param = scanPars)

# Auxilliary function to get error rate
GetConsens <- function(SMat) {
  apply(SMat, 1, FUN = function(x) {
    NucCount <- table(x)
    names(NucCount)[which.max(NucCount)]
  })
}

ErrorRegion <- function(blnDiff, WWidth = 10, ErrorThresh = 0.25){
  SLen   <- length(blnDiff)
  WStart <- seq(1, SLen, WWidth)
  WEnd   <- c(WStart[-1], SLen)
  AvDiff <- rep(NA, SLen)
  for (i in 1:length(WStart)){
    AvDiff[WStart[i]:WEnd[i]] <- mean(blnDiff[WStart[i]:WEnd[i]])
  }
  AvDiff > ErrorThresh
}

GetErrorRate <- function(RL, GR, RefSeqV, MinCoverPerZMW = 5){
  
  # Get ZMW id from read ID
  ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
    strsplit(RL$qname[j], "/")[[1]][2]
  })
  
  # Count ZMWIDs
  ZMWIDcount <- table(ZMW_IDs)
  ZMWIDcount <- ZMWIDcount[ZMWIDcount >= 5]
  if (length(ZMWIDcount) >= 2) {
    CountOrder <- order(ZMWIDcount, decreasing = T)
    ZMW2Use    <- names(ZMWIDcount)[CountOrder[1:2]]
    blnZMW     <- ZMW_IDs %in% ZMW2Use
    RLLocal    <- lapply(RL, function(y) y[blnZMW])
    
    # Get a matrix of sequences
    SeqMat <- SeqMatFromReads(RLLocal, GR)
    
    # Index for first and second ZMW
    idxZMW1 <- which(ZMW_IDs[blnZMW] == ZMW2Use[1])
    idxZMW2 <- which(ZMW_IDs[blnZMW] == ZMW2Use[2])
    
    # Get within ZMW error rate
    blnDiffWithinZMW1 <- SeqMat[,idxZMW1[1]] != SeqMat[,idxZMW1[2]]
    blnDiffWithinZMW2 <- SeqMat[,idxZMW2[1]] != SeqMat[,idxZMW2[2]]
    
    # Get the consensus sequences
    Consens1 <- GetConsens(SeqMat[,idxZMW1])
    Consens2 <- GetConsens(SeqMat[,idxZMW2])
    
    # Get consensus per ZMW ID
    blnNoStar1 <- Consens1 != "*"
    blnNoStar2 <- Consens2 != "*" 
    blnNoStar  <- blnNoStar1 & blnNoStar2
    blnDiff    <- Consens1 != Consens2
    blnMinL1   <- sum(blnNoStar1) > 500
    blnMinL2   <- sum(blnNoStar2) > 500
    blnMinL    <- sum(blnNoStar)  > 500
    blnDiffRef1 <- Consens1 != RefSeqV
    blnDiffRef2 <- Consens2 != RefSeqV
    blnError1 <- ErrorRegion(blnDiffRef1)
    blnError2 <- ErrorRegion(blnDiffRef2)
    blnError  <- blnError1 & blnError2
    
    # Calculate error rate
    c(sum(blnDiff & blnMinL & blnError & blnNoStar) / sum(blnMinL & blnError & blnNoStar),
      sum(blnDiffWithinZMW1 & blnMinL1 & blnError1 & blnNoStar1) / 
        sum(blnMinL1 & blnError1 & blnNoStar1),
      sum(blnDiffWithinZMW2 & blnMinL2 & blnError2 & blnNoStar2) / 
        sum(blnMinL2 & blnError2 & blnNoStar2),
      sum(blnDiffWithinZMW1 & blnMinL1 & (!blnError1) & blnNoStar1) / 
        sum(blnMinL1 & (!blnError1) & blnNoStar1),
      sum(blnDiffWithinZMW2 & blnMinL2 & (!blnError2) & blnNoStar2) / 
        sum(blnMinL2 & (!blnError2) & blnNoStar2))
  } else {rep(NA, 5)}
}

# Function to remove NAs
RemoveNA <- function(RL){
  blnPosNotNA <- !is.na(RL$pos)
  lapply(RL, function(x) x[blnPosNotNA])
}
cat("\n**********    Calculating error rate    *************\n")

# Initialize vectors of error rate
ErrorBetwZMW_highErr  <- c()
ErrorWithZMW_highErr  <- c()
ErrorWithZMW_Normal   <- c()
idxGRNormalBWA        <- c()

# Loop through ranges and get error rates 
for (i in 1:length(GRUnion_withSNP)){
  cat("Processing range", i, "of",  length(GRUnion_withSNP), "\n")
  # Get genomic range and reference sequence
  GR <- GRUnion_withSNP[i]
  GStart <- start(GR)
  GEnd   <- end(GR)
  GWidth <- width(GR)
  RefSeq <- RefSeqs[i]
  RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
#if (i == 98) browser()  
  # Get HiFi reads and remove NAs

  # Get Normal reads aligned with bwaand remove NAs
  RL  <- RemoveNA(ReadList_Normal_BWA[[i]])
  
  # Get ZMW id from read ID
  if (length(RL$pos) > 0){
    ErrRates  <- GetErrorRate(RL, GR, RefSeqV)
    ErrorBetwZMW_highErr <- c(ErrorBetwZMW_highErr, ErrRates[1])
    ErrorWithZMW_highErr <- c(ErrorWithZMW_highErr, list(ErrRates[2:3]))
    ErrorWithZMW_Normal  <- c(ErrorWithZMW_Normal, list(ErrRates[4:5]))
    idxGRNormalBWA    <- c(idxGRNormalBWA, i)
  }
}


# Saving data
cat("\n***** Saving data to", PathOutput, "  *****\n")
save(list = c("ErrorBetwZMW_highErr", "ErrorWithZMW_highErr", "ErrorWithZMW_Normal", 
  "idxGRNormalBWA", "GRUnion_withSNP"), 
    file = PathOutput)

