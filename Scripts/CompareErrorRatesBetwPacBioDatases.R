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
PathBam_HiFi       <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked_ngmlr.sorted.bam"
PathBam_HiFi_BWA   <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked.sorted.bam'
PathBam_Normal     <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.nglmr.sorted.bam"
PathBam_Normal_BWA <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
PathOutputNew      <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/ErrorComparer.RData'
PathPlot           <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/ErrorRateComparison.pdf"

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
ReadList_HiFi   <- scanBam(PathBam_HiFi, param = scanPars)
ReadList_HiFi_BWA   <- scanBam(PathBam_HiFi_BWA, param = scanPars)
ReadList_Normal <- scanBam(PathBam_Normal, param = scanPars)
ReadList_Normal_BWA <- scanBam(PathBam_Normal_BWA, param = scanPars)

# Auxilliary function to get error rate
GetConsens <- function(SMat) {
  apply(SMat, 1, FUN = function(x) {
    NucCount <- table(x)
    names(NucCount)[which.max(NucCount)]
  })
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
    
    # Get the consensus sequences
    Consens1 <- GetConsens(SeqMat[,ZMW_IDs[blnZMW] == ZMW2Use[1]])
    Consens2 <- GetConsens(SeqMat[,ZMW_IDs[blnZMW] == ZMW2Use[2]])
    
    # Get consensus per ZMW ID
    blnNoStar1 <- Consens1 != "*"
    blnNoStar2 <- Consens2 != "*" 
    blnNoStar  <- blnNoStar1 & blnNoStar2
    blnDiff   <- Consens1 != Consens2
    blnMinL1   <- sum(blnNoStar1) > 500
    blnMinL2   <- sum(blnNoStar2) > 500
    blnMinL   <- sum(blnNoStar) > 500
    blnDiffRef1 <- Consens1 != RefSeqV
    blnDiffRef2 <- Consens2 != RefSeqV
    
    # Calculate error rate
    c(sum(blnNoStar & blnDiff & blnMinL) / sum(blnNoStar & blnMinL),
      sum(blnNoStar1 & blnDiffRef1 & blnMinL1) / sum(blnNoStar1 & blnMinL1),
      sum(blnNoStar2 & blnDiffRef2 & blnMinL2) / sum(blnNoStar2 & blnMinL2))
  }
}

# Function to remove NAs
RemoveNA <- function(RL){
  blnPosNotNA <- !is.na(RL$pos)
  lapply(RL, function(x) x[blnPosNotNA])
}
cat("\n**********    Calculating error rate    *************\n")

# Initialize vectors of error rate
ErrorHiFi      <- c()
ErrorHiFiBWA   <- c()
ErrorNormal    <- c()
ErrorNormalBWA <- c()
ErrorRefHiFi      <- c()
ErrorRefHiFiBWA   <- c()
ErrorRefNormal    <- c()
ErrorRefNormalBWA <- c()
GRHiFi         <- c()
GRHiFiBWA      <- c()
GRNormal       <- c()
GRNormalBWA    <- c()

# Loop through ranges and get error rates 
length(GRUnion_withSNP)
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
  RL  <- RemoveNA(ReadList_HiFi[[i]])
  
  # Get error rates per ZMW
  if (length(RL$pos) > 0){
    ErrRates  <- GetErrorRate(RL, GR, RefSeqV)
    ErrorHiFi <- c(ErrorHiFi, ErrRates[1])
    ErrorRefHiFi <- c(ErrorRefHiFi, ErrRates[2:3])
    GRHiFi    <- c(GRHiFi, GR)
  }
  
  # Get HiFi BWA reads and remove NAs
  RL  <- RemoveNA(ReadList_HiFi_BWA[[i]])

  # Get error rates per ZMW
  if (length(RL$pos) > 0){
    ErrRates  <- GetErrorRate(RL, GR, RefSeqV)
    ErrorHiFiBWA <- c(ErrorHiFiBWA, ErrRates[1])
    ErrorRefHiFiBWA <- c(ErrorRefHiFiBWA, ErrRates[2:3])
    GRHiFiBWA    <- c(GRHiFiBWA, GR)
  }

  # Get Normal reads and remove NAs
  RL  <- RemoveNA(ReadList_Normal[[i]])
  
  # Get ZMW id from read ID
  if (length(RL$pos) > 0){
    ErrRates  <- GetErrorRate(RL, GR, RefSeqV)
    ErrorNormal <- c(ErrorNormal, ErrRates[1])
    ErrorRefNormal <- c(ErrorRefNormal, ErrRates[2:3])
    GRNormal    <- c(GRNormal, GR)
  }
  
  # Get Normal reads aligned with bwaand remove NAs
  RL  <- RemoveNA(ReadList_Normal_BWA[[i]])
  
  # Get ZMW id from read ID
  if (length(RL$pos) > 0){
    ErrRates  <- GetErrorRate(RL, GR, RefSeqV)
    ErrorNormalBWA <- c(ErrorNormalBWA, ErrRates[1])
    ErrorRefNormalBWA <- c(ErrorRefNormalBWA, ErrRates[2:3])
    GRNormalBWA    <- c(GRNormalBWA, GR)
  }
}
ErrorHiFi
ErrorHiFiBWA
ErrorNormal
ErrorNormalBWA
ErrorRefHiFi      
ErrorRefHiFiBWA   
ErrorRefNormal    
ErrorRefNormalBWA 
mean(ErrorHiFi, na.rm = T)
mean(ErrorHiFiBWA, na.rm = T)
mean(ErrorNormal, na.rm = T)
mean(ErrorNormalBWA, na.rm = T)
mean(ErrorRefNormal, na.rm = T)
mean(ErrorRefNormalBWA, na.rm = T)

# Plot distribution of error rates
cat("\n***** Writing plot to", PathPlot, "  *****\n")
CountBreaks        <- seq(0, 1, 0.05)
HistErrorHiFi         <- hist(ErrorHiFi,       plot = F, breaks = CountBreaks)
HistErrorHiFiBWA      <- hist(ErrorHiFiBWA,    plot = F, breaks = CountBreaks)
HistErrorNormal       <- hist(ErrorNormal,     plot = F, breaks = CountBreaks)
HistErrorNormalBWA    <- hist(ErrorNormalBWA,  plot = F, breaks = CountBreaks)
HistErrorRefHiFi      <- hist(ErrorRefHiFi,    plot = F, breaks = CountBreaks)   
HistErrorRefHiFiBWA   <- hist(ErrorRefHiFiBWA, plot = F, breaks = CountBreaks)  
HistErrorRefNormal    <- hist(ErrorRefNormal,  plot = F, breaks = CountBreaks)  
HistErrorRefNormalBWA <- hist(ErrorRefNormalBWA, plot = F, breaks = CountBreaks)

OffSet <- 3*10^-3
Cols <- rainbow(8)
pdf(file = PathPlot)
plot(HistErrorHiFi$mids - OffSet, HistErrorHiFi$density, type = "s", col = Cols[1],
     xlim = c(0, 1), xlab = "Error rate", ylab = "Frequency")
lines(HistErrorHiFiBWA$mids, HistErrorHiFiBWA$density, type = "s", col =  Cols[2])
lines(HistErrorNormal$mids, HistErrorNormal$density, type = "s", col =  Cols[3])
lines(HistErrorNormalBWA$mids, HistErrorNormalBWA$density, type = "s", col =  Cols[4])
lines(HistErrorRefHiFi$mids + OffSet, HistErrorRefHiFi$density, type = "s", col =  Cols[5])
lines(HistErrorRefHiFiBWA$mids + OffSet, HistErrorRefHiFiBWA$density, type = "s", col =  Cols[6])
lines(HistErrorRefNormal$mids + OffSet, HistErrorRefNormal$density, type = "s", col =  Cols[7])
lines(HistErrorRefNormalBWA$mids + OffSet, HistErrorRefNormalBWA$density, type = "s", col = Cols[8])
legend(x = 0.5, y = 15, legend = c("New polymerase nglmr", "New polymerase bwa", 
                                  "Old polymerase nglmr", "Old polymerase bwa",
                                  "New polymerase nglmr ref", "New polymerase bwa ref", 
                                  "Old polymerase nglmr ref", "Old polymerase bwa ref"), 
       col = Cols,
       lty = rep(1, 8), bty = "n")
dev.off()

# Saving data
cat("\n***** Saving data to", PathOutputNew, "  *****\n")
save(list = c("ErrorHiFi", "ErrorHiFiBWA", "ErrorNormal", "ErrorNormalBWA", 
 "ErrorRefHiFi", "ErrorRefHiFiBWA", "ErrorRefNormal", "ErrorRefNormalBWA",
 "GRHiFi",  "GRHiFiBWA", "GRNormalBWA"), 
    file = PathOutputNew)

