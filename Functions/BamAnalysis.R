##############################################
#
# General description:
#
#   The following function analyzes all bam and vcf files in a folder
#   

# Input:
#
#     BamFolder: character string providing path to folder that contains 
#          the fastq files to be mapped
#     L1Ranges: character string providing path to file containing the 
#          ranges of peaks intersecting with L1
#     IndexCommand: text string providing command to create an index file
#     AlignCommand: text string providing the alignment command (options can be
#          added here)
#     SamSuffix: suffix for sam files created by alignment


# Output:
#   
#     ...

##############################################

BamAnalysis <- function(BamFolder, 
                        BamSuffix,
                        L1Ranges,
                        L1Length = 6064,
                        blnGetSplitReads = F,
                        BamFileOriginal = '/home/hzudohna/NA12878-L15P_S1_L001_001.sorted.dedup.mapnonzero.bam') {
    
  cat("\n\n*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Running function BamAnalysis ...               **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")
  
  ######################################
  #                                    #
  #   Get file names and load data     #
  #                                    #
  ######################################
  
  # Load reanges of peaks intersecting with L1
  load(L1Ranges)
  
  # Get all paths to bam files in the folder
  BamFileNames <- list.files(BamFolder, pattern = BamSuffix, full.names = T)
  Bam.Pattern <- grep(".bam.", BamFileNames)
  if (length(Bam.Pattern) > 0){
    BamFileNames <- BamFileNames[-grep(".bam.", BamFileNames)]
  }
  if(length(BamFileNames) != length(idxSuspectL1Ranges)){
    stop("Number of bam files in folder does not match number of peaks!\n")
  }
  

  # Get indices from bam and vcf files
  idxFromBam <- sapply(BamFileNames, function(x){
    PathSplit <- strsplit(x, "/")[[1]]
    Fname <- PathSplit[length(PathSplit)]
    as.numeric(strsplit(Fname, "_")[[1]][2])
  })

  # Reorder bam and vcf file names so that they match idxSuspectL1Ranges
  idxMatchBam <- match(idxSuspectL1Ranges, idxFromBam)
  if (any(is.na(idxMatchBam))){
    stop("Bam file names don't match idxSuspectL1Ranges!\n")
  }
  BamFileNames <- BamFileNames[idxMatchBam]

  # Loop through file names and read in bam files of reads mapped to L1
  ScannedL1Ranges <- lapply(BamFileNames, function(x) scanBam(x))
  
  ######################################
  #                                    #
  #   Calculate basic quantities       #
  #                                    #
  ######################################
  
  # Count the number of reads mapped
  NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
    sum(!is.na(x[[1]]$pos))
  })
  idxNonZero <- which(NrMapped2L1 > 0)
  
  # Create a coverage matrix that gives the coverage per base
  CoverMat <- matrix(0, nrow = length(BamFileNames), ncol = L1Length)
  rownames(CoverMat) <- BamFileNames
  R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
                ranges = IRanges(start = 1, end = L1Length))
  for (i in idxNonZero){
    Reads <- extractReads(BamFileNames[i], R1)
    Cov <- coverage(Reads)
    CoverMat[i,] <- as.vector(Cov[[1]])
  }
  
  # Get means and 95% quantiles 
  QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
  idxFw <- 1:ncol(CoverMat)
  idxRv <- ncol(CoverMat):1
  
  ######################################
  #                                    #
  #   Summary stastistics per peak     #
  #                                    #
  ######################################
  
  # Create a data frame that keeps track of all the peak following summaries:
  # FirstPos:   First position with nonzero coverage
  # LastPos:    last position with nonzero coverage
  # MaxCovL1:   maximum coverage on L1
  # MinCovL1:   minimum coverage on L1
  # MeanCovL1:  mean coverage on L1
  # MeanCovL1NonZero:  mean coverage on L1 per position with non-zero coverage
  # MeanCovL1Range:   mean coverage on L1 per position within range with 
  #                   non-zero coverage
  # MeanCovL1:  mean coverage on L1
  # MaxCovRef:  maximum coverage on reference genome
  # WidthRef:   peak width on reference genome
  # wMismatch: Proportion different from L1 weighted by coverage
  # L1Class:    L1 classification
  EmptyData <- rep(NA, length(BamFileNames))
  PeakSummary <- data.frame(BamFile = BamFileNames,
                            FirstPos = EmptyData, 
                            LastPos = EmptyData,
                            L1Width = EmptyData,
                            MinCovL1 = apply(CoverMat, 1, min),
                            MaxCovL1 = apply(CoverMat, 1, max),
                            MeanCovL1 = rowMeans(CoverMat),
                            MeanCovL1NonZero = EmptyData,
                            MeanCovL1Range = EmptyData,
                            MaxCovRef = maxCover[idxSuspectL1Ranges],
                            WidthRef  = width(SuspectL1Ranges),
                            wMismatch = EmptyData,
                            L1Class = EmptyData,
                            idxSuspectL1Ranges = idxSuspectL1Ranges,
                            startRef = start(IslGRanges_reduced[idxSuspectL1Ranges]),
                            endRef = end(IslGRanges_reduced[idxSuspectL1Ranges])
  )
  for (i in 1:nrow(PeakSummary)){
    Cov        <- CoverMat[i,]
    NonzeroCov <- which(Cov > 0)
    if (sum(Cov) > 0){
      PeakSummary$FirstPos[i] <- min(NonzeroCov)
      PeakSummary$LastPos[i]  <- max(NonzeroCov)
      PeakSummary$MeanCovL1NonZero[i] <- mean(Cov[NonzeroCov])
      PeakSummary$MeanCovL1Range[i] <- mean(Cov[min(NonzeroCov):max(NonzeroCov)])
      PeakSummary$L1Width[i]  <- max(NonzeroCov) - min(NonzeroCov)
    }
  }
  
  if (blnGetSplitReads){
    # Loop through bamfiles and determine indices of reads that
    # were mapped to L1 and the reference genome
    SplitReadList <- lapply(BamFileNames, function(BamFile){
      
      # Get reads mapped to reference genome
      idxRange      <- as.numeric(strsplit(BamFile, "_")[[1]][3])
      IslandRange   <- IslGRanges_reduced[idxRange]
      param         <- ScanBamParam(which = IslandRange, what = scanBamWhat())
      OriginalReads <- scanBam(BamFileOriginal, param = param)
      
      # Get reads for L1 ranges
      ScannedReads <- scanBam(BamFile)
      
      # Match IDs
      IDMAtch <- match(ScannedReads[[1]]$qname, OriginalReads[[1]]$qname)
      OrigninalStrand <- OriginalReads[[1]]$strand[IDMAtch]
      
      # Get split reads
      idxSplitRead <- (OriginalReads[[1]]$seq[IDMAtch] == ScannedReads[[1]]$seq) &
        !is.na(ScannedReads[[1]]$pos)
      if (any(idxSplitRead)){
        Output <- list(Seqs = ScannedReads[[1]]$seq[idxSplitRead],
                       CigarOrig = OriginalReads[[1]]$cigar[IDMAtch][idxSplitRead],
                       CigarL1 = ScannedReads[[1]]$cigar[idxSplitRead],
                       PosOrig = OriginalReads[[1]]$pos[IDMAtch][idxSplitRead],
                       PosL1 = ScannedReads[[1]]$pos[idxSplitRead])
      } else {
        Output <- NULL
      }
      Output
    })
    
  } else {
    SplitReadList <- NULL
  }
  list(CoverMat = CoverMat, QuantileMat = QuantileMat, PeakSummary = PeakSummary,
       PeakRanges = IslGRanges_reduced[idxSuspectL1Ranges],
       SplitReadList = SplitReadList)
}


