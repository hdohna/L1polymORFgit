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
# Paths (local computer)
StartPath <- 'D:/L1polymORFgit/Scripts/_Start_L1polymORF.R'
DataPath  <- "D:/L1polymORF/Data/ErrorComparer.RData"
PathBam_Normal_BWA <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19masked_highErrorchr8.bam"
PathBam_Normal_BWA <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19masked_highErrorchr18.bam"
PathBam_Normal_BWA <- "D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19masked_highErrorchr16.bam"

# Source start script
source(StartPath)

###############################################
#                                             #
#    Load ranges where reads were mapped      #
#                                             #
###############################################

cat("\n**********    Loading data    *************\n")

# Get example reads from a low fidelity region
ReadList_Normal_BWA <- scanBam(PathBam_Normal_BWA)

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

DNAStFromSeqMat <- function(SeqMat, seqSubset, posSubset){
  Seq <- SeqMat[posSubset, seqSubset]
  Seq <- Seq[!Seq %in% c("-", "*")]
  DNAString(paste(Seq, collapse = ""))
}
DNAStFromSeq <- function(Seq){
  Seq <- Seq[!Seq %in% c("-", "*")]
  DNAString(paste(Seq, collapse = ""))
}

RL <- ReadList_Normal_BWA[[1]]
GRs <- ReadList2GRanges(RL)
GR <- union(GRs, GRs)
RefSeq  <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GR)
RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]

#GetErrorRate <- function(RL, GR, RefSeqV, MinCoverPerZMW = 5){
  
  # Get ZMW id from read ID
  ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
    strsplit(RL$qname[j], "/")[[1]][2]
  })
  
  # Count ZMWIDs
  ZMWIDcount <- table(ZMW_IDs)
  ZMWIDcount <- ZMWIDcount[ZMWIDcount >= 5]
#  if (length(ZMWIDcount) >= 2) {
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
    
    # Get sequences in error regions
    Seq1DNA_Err <- DNAStFromSeqMat(SeqMat, seqSubset = idxZMW1[1],
                               posSubset = blnError & blnNoStar)
    Seq2DNA_Err <- DNAStFromSeqMat(SeqMat, seqSubset = idxZMW2[3],
                                   posSubset = blnError & blnNoStar)
    Seq1DNA_NotErr <- DNAStFromSeqMat(SeqMat, seqSubset = idxZMW1[1],
                                   posSubset = (!blnError) & blnNoStar)
    Seq2DNA_NotErr <- DNAStFromSeqMat(SeqMat, seqSubset = idxZMW1[3],
                                   posSubset = (!blnError) & blnNoStar)
    PWA_Err <- pairwiseAlignment(Seq1DNA_Err, Seq2DNA_Err)
    PDiff_Err <- (length(PWA_Err@pattern@mismatch[[1]]) +
        sum(PWA_Err@pattern@indel[[1]]@width) + sum(PWA_Err@subject@indel[[1]]@width))/
        width(PWA_Err@pattern)
    PWA_NotErr <- pairwiseAlignment(Seq1DNA_NotErr, Seq2DNA_NotErr)
    PDiff_NotErr <- (length(PWA_NotErr@pattern@mismatch[[1]]) +
                    sum(PWA_NotErr@pattern@indel[[1]]@width) + sum(PWA_NotErr@subject@indel[[1]]@width))/
      width(PWA_NotErr@pattern)
    PDiff_Err
    PDiff_NotErr
    # Calculate error rate
 #}

