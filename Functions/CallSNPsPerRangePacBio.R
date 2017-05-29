##############################################
#
# General description:
#
#   The following function takes a read list (as created by the function 
#   scanBam), a genomic range and calls all SNPs within that range

# Input:
#
#     RL: list with read info (as created by the function scanBam)
#     GR: Genomic range within which to call SNPs
#     MinProp: minimum proportion of reads that should contain SNP

# Output:
#   
#    SNPpos: SNP positions

# Comments:
#   
#    Requires function SeqFromCigar and ReadList2GRanges

##############################################

CallSNPsPerRangePacBio <- function(RL, GR, RefGenome = BSgenome.Hsapiens.UCSC.hg19,
                                   MinProp = 0.6){
  
  # Set parameters of current range
  GStart  <- start(GR)
  GEnd    <- end(GR)
  GWidth  <- width(GR)
  RefSeq  <- getSeq(RefGenome, GR)
  RefSeqV <- strsplit(as.character(RefSeq), "")[[1]]
  Chr     <- RL$rname[1]
  
  # Subset read list to get only the reads intersecting with current range
  # (Maybe move outside function)
  ReadGR     <- ReadList2GRanges(RL)
  blnOverlap <- overlapsAny(ReadGR, GR)
  RL         <- lapply(RL, function(x) x[blnOverlap])
  
  # Get a list of all aligned sequences
  SeqList <- lapply(1:length(RL$pos), function(j) {
    SeqFromCigar(RL$cigar[j], RL$seq[j])
  })

  # Get all positions that differ from reference
  DiffPosList <- lapply(1:length(RL$pos), function(j) {
      SeqV     <- SeqList[[j]]
      SeqStart <- max(1, GStart - RL$pos[j] + 1)
      SeqEnd   <- min(GEnd - RL$pos[j] + 1, length(SeqV))
      RefStart <- max(1, RL$pos[j] - GStart + 1)
      RefEnd   <- min(RL$pos[j] + length(SeqV) - GStart, GWidth)
      if(length(SeqStart:SeqEnd) != length(RefStart:RefEnd)) browser()
      RL$pos[j] + SeqStart + which(SeqV[SeqStart:SeqEnd] != RefSeqV[RefStart:RefEnd]) - 1
    })
  DiffPos <- unlist(DiffPosList)
  DiffPos <- DiffPos[duplicated(DiffPos)]
    
  # Get ZMW id from read ID
  ZMW_IDs <- sapply(1:length(RL$pos), function(j) {
      strsplit(RL$qname[j], "/")[[1]][2]
  })
  ZMW_Count <- table(ZMW_IDs)
  
  # Get positions that differ per ZMW ID
  DiffPosList_PerZMW <- lapply(ZMW_IDs, function(x){
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
         DiffPosAll = DiffPosAll,
         NrDiffPosInAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         DiffPosAll_SNP = sum(as.numeric(DiffPosInAll) %in% SNPposRange),
         NrSNPs = length(SNPpos))
}