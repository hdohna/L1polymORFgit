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
PathStart        <- '/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R'
PathSNPRanges    <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/SNPinRangesNA12878.RData'
PathBam_HiFi     <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked_ngmlr.sorted.bam"
PathBam_HiFi_BWA <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked.sorted.bam'
PathBam_Normal   <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
PathOutputNew    <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/ErrorComparer.RData'
PathPlot         <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/ErrorRateComparison.pdf"

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

# Auxilliary function to get error rate
GetErrorRate <- function(RL, RefSeqV){
  
  # Get ZMW id from read ID
  ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
    strsplit(RL$qname[j], "/")[[1]][2]
  })
  
  # Get consensus per ZMW ID
  ErrorRate <- sapply(unique(ZMW_IDs), function(x) {
    blnZMW  <- ZMW_IDs == x
    RLLocal <- lapply(RL, function(y) y[blnZMW])
    ConsSeq <- ConsensusFromReads(RLLocal, GR)
    Rate    <- sum(ConsSeq != RefSeqV & ConsSeq != "*") / sum(ConsSeq != "*")
    c(Rate, NA)[1 + ((sum(blnZMW) < 3) | (sum(ConsSeq != "*") < 500))]
  })
  ErrorRate[!is.na(ErrorRate)]
}

cat("\n**********    Calculating error rate    *************\n")

# Initialize vectors of error rate
ErrorHiFi     <- c()
ErrorHiFiBWA  <- c()
ErrorNormal   <- c()
GRHiFi        <- c()
GRHiFiBWA     <- c()
GRNormal      <- c()

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
  
  # Get HiFi reads and remove NAs
  RL          <- ReadList_HiFi[[i]]
  blnPosNotNA <- !is.na(RL$pos)
  RL          <- lapply(RL, function(x) x[blnPosNotNA])
  
  # Get error rates per ZMW
  if (length(RL$pos) > 0){
    ErrorHiFi <- c(ErrorHiFi, GetErrorRate(RL, RefSeqV))
    GRHiFi    <- c(GRHiFi, GR)
  }
  
  # Get HiFi BWA reads and remove NAs
  RL          <- ReadList_HiFi_BWA[[i]]
  blnPosNotNA <- !is.na(RL$pos)
  RL          <- lapply(RL, function(x) x[blnPosNotNA])
  
  # Get error rates per ZMW
  if (length(RL$pos) > 0){
    ErrorHiFiBWA <- c(ErrorHiFiBWA, GetErrorRate(RL, RefSeqV))
    GRHiFiBWA    <- c(GRHiFiBWA, GR)
  }

    # Get Normal reads and remove NAs
  RL          <- ReadList_Normal[[i]]
  blnPosNotNA <- !is.na(RL$pos)
  RL          <- lapply(RL, function(x) x[blnPosNotNA])
  
  # Get ZMW id from read ID
  if (length(RL$pos) > 0){
    ErrorNormal <- c(ErrorNormal, GetErrorRate(RL, RefSeqV))
    GRNormal    <- c(GRNormal, GR)
  }
  
}
ErrorHiFi
ErrorHiFiBWA
ErrorNormal
mean(ErrorHiFi)
mean(ErrorHiFiBWA)
mean(ErrorNormal)

# Plot distribution of error rates
cat("\n***** Writing plot to", PathPlot, "  *****\n")
CountBreaks     <- seq(0, 1, 0.05)
HistErrorHiFi   <- hist(ErrorHiFi,   plot = F, breaks = CountBreaks)
HistErrorHiFiBWA   <- hist(ErrorHiFiBWA,   plot = F, breaks = CountBreaks)
HistErrorNormal <- hist(ErrorNormal, plot = F, breaks = CountBreaks)
OffSet <- 3*10^-3
pdf(file = PathPlot)
plot(HistErrorHiFi$mids - OffSet, HistErrorHiFi$density, type = "s", col = "blue",
     xlim = c(0, 1), xlab = "Error rate", ylab = "Frequency")
lines(HistErrorHiFiBWA$mids, HistErrorHiFiBWA$density, type = "s", col = "green")
lines(HistErrorNormal$mids + OffSet, HistErrorNormal$density, type = "s", col = "red")
legend(x = 0.6, y = 9, legend = c("New polymerase nglmr", "New polymerase bwa", 
                                  "Old polymerase bwa"), col = c("blue", "green", "red"),
       lty = c(1, 1, 1), bty = "n")
dev.off()

cat("\n***** Saving data to", PathOutputNew, "  *****\n")
save(list = c("ErrorHiFi", "GRHiFi", "ErrorHiFiBWA", "GRHiFiBWA", "ErrorNormal", "GRNormal"), 
     file = PathOutputNew)

