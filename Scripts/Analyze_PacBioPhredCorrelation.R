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
PathStart      <- '/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R'
PathSNPRanges  <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/SNPinRangesNA12878.RData'
PathOutputNew     <- '/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PhredCor.RData'

# Source start script
source(PathStart)

###############################################
#                                             #
#    Load ranges where reads were mapped      #
#                                             #
###############################################

# Load data with ranges wehere reads were mapped and SNPs
load(PathSNPRanges)

#ReadList_HiFi   <- scanBam(PathBam_HiFi)

# Get ranges with minimum coverage
MinCoverRanges_HiFi <- GRanges()
length(ReadList_HiFi)

PhredCor_sameZMW <- c()
PhredCor_diffZMW <- c()
PhredCor_Random  <- c()
LowPhredMotif    <- c()
LowPhredMotifRandom    <- c()

#for (i in 1:length(ReadList_HiFi)){
for (i in 1:length(ReadList_HiFi)){
  print(i)
  RL <- ReadList_HiFi[[i]]
  if (length(RL$qname) > 0){
    blnNaPos <- is.na(RL$pos)
    RL <- lapply(RL, function(x) x[!blnNaPos])
    RGR      <- ReadList2GRanges(RL)
    RCover   <- coverage(RGR)[RL$rname[1]][[1]]
    RCoverIR <- slice(RCover, lower =  5, rangesOnly = T)
    RCoverIR <- RCoverIR[which.max(width(RCoverIR))]
    Rseqname <- rep(RL$rname[1], length(RCoverIR))
    RCoverGR <- GRanges(Rseqname, ranges = RCoverIR)
    
    if (length(RCoverGR) > 0) {
      # Get a list of all aligned phred scores
      StartAll   <- start(RCoverIR)
      EndAll     <- end(RCoverIR)
      blnInRange <- overlapsAny(RGR, RCoverGR,  minoverlap = width(RCoverGR))
      RL <- lapply(RL, function(x) x[blnInRange])
      
      if (length(RL$pos) > 0) {
        # Get a matrix of phred scores
        PhredMat <- sapply(seq_along(RL$pos), function(j) {
          PhredV     <- PhredFromCigar(RL$cigar[j], RL$qual[j])
          PhredStart <- StartAll - RL$pos[j] + 1
          PhredEnd   <- EndAll - RL$pos[j] + 1
          PhredV[PhredStart:PhredEnd]
        })
        
        # Calculate a correlation of phred scores
        PhredCor <- cor(PhredMat, use = "pairwise.complete.obs")
        diag(PhredCor) <- NA
        
        # Calculate a correlation between random phred scores
        PhredMatRandom <- sapply(1:ncol(PhredMat), function(x){
          idxRow <- sample(1:nrow(PhredMat), nrow(PhredMat))
          PhredMat[idxRow, x]
        })
        PhredCorRandom <- cor(PhredMatRandom, use = "pairwise.complete.obs")
        diag(PhredCorRandom) <- NA
        
        # Get ZMW id from read ID
        ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
          strsplit(RL$qname[j], "/")[[1]][2]
        })
        blnSameZMW_ID <- outer(ZMW_IDs, ZMW_IDs, function(x, y) x == y)
        
        # Get the position of lowest mean phred score
        PhredMean <- rowMeans(PhredMat)
        MinPos    <- which.min(PhredMean)
        MotifRange <- IRanges(start = StartAll + max(0, MinPos - 2), 
                              end = min(EndAll, StartAll + MinPos + 2))
        MotifRange <- GRanges(RL$rname[1], ranges = MotifRange)
        
        # Get the position of lowest mean phred score
        PhredMeanRandom <- rowMeans(PhredMatRandom)
        MinPosRandom    <- which.min(PhredMeanRandom)
        MotifRangeRandom <- IRanges(start = StartAll + max(0, MinPosRandom - 2), 
                              end = min(EndAll, StartAll + MinPosRandom + 2))
        MotifRangeRandom <- GRanges(RL$rname[1], ranges = MotifRangeRandom)
        
        # Get a matrix of phred scores
        PhredMat <- sapply(1:length(RL$pos), function(j) {
          PhredV     <- PhredFromCigar(RL$cigar[j], RL$qual[j])
          PhredStart <- StartAll - RL$pos[j] + 1
          PhredEnd   <- EndAll - RL$pos[j] + 1
          PhredV[PhredStart:PhredEnd]
        })
        
        # Update quantities
        PhredCor_sameZMW <- c(PhredCor_sameZMW,  mean(PhredCor[blnSameZMW_ID], na.rm = T))
        PhredCor_diffZMW <- c(PhredCor_diffZMW,  mean(PhredCor[!blnSameZMW_ID], na.rm = T))
        PhredCor_Random  <- c(PhredCor_Random,  mean(PhredCorRandom, na.rm = T))
        LowPhredMotif <- c(LowPhredMotif, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, MotifRange)[[1]]))
        LowPhredMotifRandom <- c(LowPhredMotifRandom, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, MotifRangeRandom)[[1]]))
        
      }
       
    }   
  }
}

# Save quantities
PhredCor_sameZMW 
PhredCor_diffZMW 
PhredCor_Random  
LowPhredMotif    
unlist(LowPhredMotif)
mean(c(PhredCor_sameZMW, PhredCor_diffZMW), na.rm = T)
mean(PhredCor_Random, na.rm = T)
mean(PhredCor_sameZMW, na.rm = T)
mean(PhredCor_diffZMW, na.rm = T)

LowPhredMotifCount <- table(LowPhredMotif)
LowPhredMotifCount[order(LowPhredMotifCount)]

LowPhredMotifCountRandom <- table(LowPhredMotifRandom)
LowPhredMotifCountRandom[order(LowPhredMotifCountRandom)]
NameMatch <- match(names(LowPhredMotifCount), names(LowPhredMotifCountRandom))
CountRatio <- LowPhredMotifCount[!is.na(NameMatch)] / 
  LowPhredMotifCountRandom[NameMatch[!is.na(NameMatch)]]
sort(log(CountRatio))

# Save objects indicating phred correlation
save(list = c("PhredCor_sameZMW", "PhredCor_diffZMW", "PhredCor_Random", "LowPhredMotif", 
              "LowPhredMotifRandom"), file = PathOutputNew)
