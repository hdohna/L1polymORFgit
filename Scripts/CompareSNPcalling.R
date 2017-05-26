# This script compares SNP calling between two datasets

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

# Paths (local computer)
PathStart             <- 'D:/L1polymORFgit/Scripts/_Start_L1polymORF.R'
PathReadRanges_HiFi   <- "D:/L1polymORF/Data/PacBioHiFiSubreads_L1Ranges.RData"
PathReadRanges_Normal <- "D:/L1polymORF/Data/BZ_L1Ranges.RData"
PathSNP_NA12878       <- "D:/L1polymORF/Data/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf"
PathBam_HiFi          <- ""
PathBam_Normal        <- ""

# Paths (SCG4)
PathStart             <- '/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R'
PathChrLen            <- '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
PathReadRanges_HiFi   <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioHiFiSubreads_L1Ranges.RData"
PathReadRanges_Normal <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_L1Ranges.RData"
PathSNP_NA12878       <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf"
PathBam_HiFi          <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked_ngmlr.sorted.bam"
PathBam_Normal        <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"

# Source start script
source(PathStart)

###############################################
#                                             #
#    Load ranges where reads were mapped      #
#                                             #
###############################################

# Load GRanges from HiFi polymerase data
load(PathReadRanges_HiFi)
IslGRanges_reduced_HiFi <- IslGRanges_reduced

# Load GRanges from HiFi polymerase data
load(PathReadRanges_Normal)
IslGRanges_reduced_Normal <- IslGRanges_reduced

# Get intersection of both sets of genomic range
GRUnion <- union(IslGRanges_reduced_HiFi, IslGRanges_reduced_Normal)

load(PathChrLen)

###############################################
#                                             #
#    Get SNPs overlapping with above ranges   #
#                                             #
###############################################

# Get SNPs of NA12878 overlapping with genomic ranges
VCFLines <- readLines(PathSNP_NA12878, n = 200)

SNP_GR       <- GRanges()
NrLines2Read <- 10^5
LinesSkip    <- max(grep("##", VCFLines)) + 1
NrLinesRead  <- NrLines2Read
while(NrLinesRead >= NrLines2Read){
  SNPs <- read.delim(file = PathSNP_NA12878, skip = LinesSkip, nrows = NrLines2Read,
                    header = F)
  SNPs <- SNPs[,c(1:2, 2, 4, 5)]
  colnames(SNPs)  <- c("chromosome", "start", "end", "REF", "ALT")
  SNPs$chromosome <- paste("chr", SNPs$chromosome, sep = "")
  SNP_GR_local    <- makeGRangesFromDataFrame(SNPs, keep.extra.columns = T)
  SNP_GR_local    <- subsetByOverlaps(SNP_GR_local, GRUnion)
  SNP_GR          <- c(SNP_GR, SNP_GR_local)
  NrLinesRead     <- nrow(SNPs)
  LinesSkip       <- LinesSkip + NrLinesRead
  cat("Processed", LinesSkip, "lines\n")
}
rm(list = "SNPs")
gc()

###############################################
#                                             #
#    Get SNPs overlapping with above ranges   #
#                                             #
###############################################

# Subset GRUnion to get the ones containing SNPs 
GRUnion_withSNP <- subsetByOverlaps(GRUnion, SNP_GR)

# Get reference sequence for each range
RefSeqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRUnion_withSNP)

# Get Reads for 
scanPars        <- ScanBamParam(which = GRUnion_withSNP, what = scanBamWhat())
ReadList_HiFi   <- scanBam(PathBam_HiFi, param = scanPars)
ReadList_Normal <- scanBam(PathBam_Normal, param = scanPars)

# Get ranges with minimum coverage
MinCoverRanges_HiFi <- GRanges()
for (i in 1:length(ReadList_HiFi)){
  print(i)
  RL <- ReadList_HiFi[[i]]
  if (length(RL$qname) > 0){
    blnNaPos <- is.na(RL$pos)
    RL <- lapply(RL, function(x) x[!blnNaPos])
    MCR <- MinCoverRangesPacBio(RL, ChrLen = ChromLengthsHg19[as.character(RL$rname[1])])
    MinCoverRanges_HiFi <- c(MinCoverRanges_HiFi, MCR)
    
  }
}





# Get indicator for regions with enough ZMW ids to call SNPs
i <- 4
blnEnoughZMW_HiFi <- sapply(ReadList_HiFi, function(RL){
  if (length(RL$pos) > 0){
    ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
      strsplit(RL$qname[j], "/")[[1]][2]
    })
    ZMW_Count <- table(ZMW_IDs)
    sum(ZMW_Count > 2) > 1
    
  } else {
    FALSE
  }
})

# PropSNP_inAll_HiFi <- sapply(which(blnEnoughZMW_HiFi), function(i){
#   RL <- ReadList_HiFi[[i]]
#   GStart <- start(GRUnion_withSNP)[i]
#   GEnd   <- end(GRUnion_withSNP)[i]
#   GWidth <- width(GRUnion_withSNP)[i]
#   RefSeq <- RefSeqs[i]
#   RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
#   
#   # Get a list of all 
#   SeqList <- lapply(1:length(RL$pos), function(j) {
#     SeqFromCigar(RL$cigar[j], RL$seq[j])
#   })
#   StartAll <- max(c(GStart, RL$pos))
#   EndAll   <- min(c(GEnd, RL$pos + sapply(SeqList, length) - 1))
#   RefStart <- StartAll - GStart + 1
#   RefEnd   <- GWidth - (GEnd - EndAll)
#   
#   if (RefEnd > RefStart){
#     DiffPosList <- lapply(1:length(RL$pos), function(j) {
#       SeqV     <- SeqList[[j]]
#       SeqStart <- StartAll - RL$pos[j] + 1
#       SeqEnd   <- EndAll - RL$pos[j] + 1
#       if(length(SeqStart:SeqEnd) != length(RefStart:RefEnd)) browser()
#       RL$pos[j] + SeqStart + which(SeqV[SeqStart:SeqEnd] != RefSeqV[RefStart:RefEnd]) - 1
#     })
#     DiffPos <- unlist(DiffPosList)
#     DiffPos <- DiffPos[duplicated(DiffPos)]
#     DiffPos
#     
#     # Get ZMW id from read ID
#     ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
#       strsplit(RL$qname[j], "/")[[1]][2]
#     })
#     ZMW_IDs
#     ZMW_Count <- table(ZMW_IDs)
#     ZMW_ID_subset <- names(ZMW_Count)[ZMW_Count > 2]
#     
#     DiffPosList_PerZMW <- lapply(ZMW_ID_subset, function(x){
#       idxSubset <- which(ZMW_IDs == x)
#       DiffPosSubset <- DiffPosList[idxSubset]
#       DiffPosCount  <- table(unlist(DiffPosSubset))
#       names(DiffPosCount)[(DiffPosCount / length(idxSubset)) > MinPropCall ]
#     })
#     DiffPos_PerZMW <- unlist(DiffPosList_PerZMW)
#     sum(table(DiffPos_PerZMW) == length(DiffPosList_PerZMW)) / length(unique(DiffPos_PerZMW))
#     
#   } else {
#     NULL
#   }
# })
# 
# mean(unlist(PropSNP_inAll_HiFi), na.rm = T)

DiffPosList_HiFi <- lapply(which(blnEnoughZMW_HiFi), function(i){
  RL <- ReadList_HiFi[[i]]
  GStart <- start(GRUnion_withSNP)[i]
  GEnd   <- end(GRUnion_withSNP)[i]
  GWidth <- width(GRUnion_withSNP)[i]
  RefSeq <- RefSeqs[i]
  RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
  Chr <- RL$rname[1]
  
  # Get a list of all aligned sequences
  SeqList <- lapply(1:length(RL$pos), function(j) {
    SeqFromCigar(RL$cigar[j], RL$seq[j])
  })
  StartAll <- max(c(GStart, RL$pos))
  EndAll   <- min(c(GEnd, RL$pos + sapply(SeqList, length) - 1))
  RefStart <- StartAll - GStart + 1
  RefEnd   <- GWidth - (GEnd - EndAll)
  
  
  if (RefEnd > RefStart){
    LocalGR <- GRanges(Chr, IRanges(start = StartAll, end = EndAll))
    SNPsLocal <- subsetByOverlaps(SNP_GR, LocalGR)
    SNPpos    <- start(SNPsLocal)
    SNPposRange <- unlist(lapply((-3):3, function(x) SNPpos + x))
    DiffPosList <- lapply(1:length(RL$pos), function(j) {
      SeqV     <- SeqList[[j]]
      SeqStart <- StartAll - RL$pos[j] + 1
      SeqEnd   <- EndAll - RL$pos[j] + 1
      if(length(SeqStart:SeqEnd) != length(RefStart:RefEnd)) browser()
      RL$pos[j] + SeqStart + which(SeqV[SeqStart:SeqEnd] != RefSeqV[RefStart:RefEnd]) - 1
    })
    DiffPos <- unlist(DiffPosList)
    DiffPos <- DiffPos[duplicated(DiffPos)]
    DiffPos
    
    # Get ZMW id from read ID
    ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
      strsplit(RL$qname[j], "/")[[1]][2]
    })
    ZMW_IDs
    ZMW_Count <- table(ZMW_IDs)
    ZMW_ID_subset <- names(ZMW_Count)[ZMW_Count > 2]
    
    DiffPosList_PerZMW <- lapply(ZMW_ID_subset, function(x){
      idxSubset <- which(ZMW_IDs == x)
      DiffPosSubset <- DiffPosList[idxSubset]
      DiffPosCount  <- table(unlist(DiffPosSubset))
      names(DiffPosCount)[(DiffPosCount / length(idxSubset)) > MinPropCall ]
    })
    DiffPos_PerZMW <- unlist(DiffPosList_PerZMW)
    DiffPosCount_PerZMW <- table(DiffPos_PerZMW)
    DiffPosInAll <- names(DiffPosCount_PerZMW)[DiffPosCount_PerZMW == length(DiffPosList_PerZMW)]
    DiffPosAll   <-  unique(DiffPos_PerZMW)
    list(DiffPosInAll = DiffPosInAll,
         DiffPosAll   =DiffPosAll,
         NrDiffPosInAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         DiffPosAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         NrSNPs = length(SNPpos))
    
  } else {
    NULL
  }
})

# Get indicator for regions with enough ZMW ids to call SNPs
blnEnoughZMW_Normal <- sapply(ReadList_Normal, function(RL){
  if (length(RL$pos) > 0){
    ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
      strsplit(RL$qname[j], "/")[[1]][2]
    })
    ZMW_Count <- table(ZMW_IDs)
    sum(ZMW_Count > 2) > 1
    
  } else {
    FALSE
  }
})

DiffPosList_Normal <- lapply(which(blnEnoughZMW_Normal), function(i){
  RL <- ReadList_Normal[[i]]
  GStart <- start(GRUnion_withSNP)[i]
  GEnd   <- end(GRUnion_withSNP)[i]
  GWidth <- width(GRUnion_withSNP)[i]
  RefSeq <- RefSeqs[i]
  RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
  Chr <- RL$rname[1]
  
  # Get a list of all 
  SeqList <- lapply(1:length(RL$pos), function(j) {
    SeqFromCigar(RL$cigar[j], RL$seq[j])
  })
  StartAll <- max(c(GStart, RL$pos))
  EndAll   <- min(c(GEnd, RL$pos + sapply(SeqList, length) - 1))
  RefStart <- StartAll - GStart + 1
  RefEnd   <- GWidth - (GEnd - EndAll)
  
  
  if (RefEnd > RefStart){
    LocalGR <- GRanges(Chr, IRanges(start = StartAll, end = EndAll))
    SNPsLocal <- subsetByOverlaps(SNP_GR, LocalGR)
    SNPpos    <- start(SNPsLocal)
    SNPposRange <- unlist(lapply((-3):3, function(x) SNPpos + x))
    DiffPosList <- lapply(1:length(RL$pos), function(j) {
      SeqV     <- SeqList[[j]]
      SeqStart <- StartAll - RL$pos[j] + 1
      SeqEnd   <- EndAll - RL$pos[j] + 1
      if(length(SeqStart:SeqEnd) != length(RefStart:RefEnd)) browser()
      RL$pos[j] + SeqStart + which(SeqV[SeqStart:SeqEnd] != RefSeqV[RefStart:RefEnd]) - 1
    })
    DiffPos <- unlist(DiffPosList)
    DiffPos <- DiffPos[duplicated(DiffPos)]
    DiffPos
    
    # Get ZMW id from read ID
    ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
      strsplit(RL$qname[j], "/")[[1]][2]
    })
    ZMW_IDs
    ZMW_Count <- table(ZMW_IDs)
    ZMW_ID_subset <- names(ZMW_Count)[ZMW_Count > 2]
    
    DiffPosList_PerZMW <- lapply(ZMW_ID_subset, function(x){
      idxSubset <- which(ZMW_IDs == x)
      DiffPosSubset <- DiffPosList[idxSubset]
      DiffPosCount  <- table(unlist(DiffPosSubset))
      names(DiffPosCount)[(DiffPosCount / length(idxSubset)) > MinPropCall ]
    })
    DiffPos_PerZMW <- unlist(DiffPosList_PerZMW)
    DiffPosCount_PerZMW <- table(DiffPos_PerZMW)
    DiffPosInAll <- names(DiffPosCount_PerZMW)[DiffPosCount_PerZMW == length(DiffPosList_PerZMW)]
    DiffPosAll   <-  unique(DiffPos_PerZMW)
    list(DiffPosInAll = DiffPosInAll,
         DiffPosAll   =DiffPosAll,
         NrDiffPosInAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         DiffPosAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         NrSNPs = length(SNPpos))
    
  } else {
    NULL
  }
})

# MeanProps
MeanProps <- c(HiFi = mean(unlist(PropSNP_inAll_HiFi), na.rm = T),
  Normal = mean(unlist(PropSNP_inAll_Normal), na.rm = T))

# Get the % of SNPs detected
NrSNPs_HiFi <- 0
NrSNPsCalled_HiFi <- 0
NrSNPsCorrectCalled_HiFi <- 0
DPList <- DiffPosList_HiFi[[1]]
for (DPList in DiffPosList_HiFi){
  if (length(DPList) > 0){
    NrSNPs_HiFi <- NrSNPs_HiFi + DPList$NrSNPs
    NrSNPsCalled_HiFi <- NrSNPsCalled_HiFi + length(DPList$DiffPosInAll)
    NrSNPsCorrectCalled_HiFi <- NrSNPsCorrectCalled_HiFi + DPList$DiffPosAll_SNP
  }
}
NrSNPsCorrectCalled_HiFi / NrSNPs_HiFi
NrSNPsCorrectCalled_HiFi / NrSNPsCalled_HiFi

NrSNPs_Normal <- 0
NrSNPsCalled_Normal <- 0
NrSNPsCorrectCalled_Normal <- 0
for (DPList in DiffPosList_Normal){
  if (length(DPList) > 0){
    NrSNPs_Normal <- NrSNPs_Normal + DPList$NrSNPs
    NrSNPsCalled_Normal <- NrSNPsCalled_Normal + length(DPList$DiffPosInAll)
    NrSNPsCorrectCalled_Normal <- NrSNPsCorrectCalled_Normal + DPList$DiffPosAll_SNP
  }
}
NrSNPsCorrectCalled_Normal / NrSNPs_Normal
NrSNPsCorrectCalled_Normal / NrSNPsCalled_Normal

