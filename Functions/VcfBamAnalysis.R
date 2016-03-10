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

VCfBamAnalysis <- function(BamFolder, 
                           BamSuffix,
                           L1Ranges,
                           L1Length = 6064,
                           CoverSummaryPlot,
                           CoverComparePlot) {
    
  cat("*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Running function VCfBamAnalysis ...            **\n")
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
  
  # Get all paths to vcf files in the folder
  VcfFileNames <- list.files(BamFolder, pattern = ".vcf", full.names = T)
  Vcf.Pattern <- grep(".vcf.", VcfFileNames)
  if (length(Bam.Pattern) > 0){
    VcfFileNames <- VcfFileNames[-grep(".vcf.", VcfFileNames)]
  }
  if(length(VcfFileNames) != length(idxSuspectL1Ranges)){
    stop("Number of vcf files in folder does not match number of peaks!\n")
  }
  
  # Get indices from bam and vcf files
  idxFromBam <- sapply(BamFileNames, function(x){
    PathSplit <- strsplit(x, "/")[[1]]
    Fname <- PathSplit[length(PathSplit)]
    as.numeric(strsplit(Fname, "_")[[1]][2])
  })
  idxFromVcf <- sapply(VcfFileNames, function(x){
    PathSplit <- strsplit(x, "/")[[1]]
    Fname <- PathSplit[length(PathSplit)]
    as.numeric(strsplit(Fname, "_")[[1]][2])
  })
  
  # Reorder bam and vcf file names so that they match idxSuspectL1Ranges
  idxMatchBam <- match(idxSuspectL1Ranges, idxFromBam)
  if (any(is.na(idxMatchBam))){
    stop("Bam file names don't match idxSuspectL1Ranges!\n")
  }
  idxMatchVcf <- match(idxSuspectL1Ranges, idxFromVcf)
  if (any(is.na(idxMatchVcf))){
    stop("Vcf file names don't match idxSuspectL1Ranges!\n")
  }
  BamFileNames <- BamFileNames[idxMatchBam]
  VcfFileNames <- VcfFileNames[idxMatchVcf]
  
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
                            VcfFile = VcfFileNames,
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
                            L1Class = EmptyData)
  for (i in 1:nrow(PeakSummary)){
    Cov        <- CoverMat[i,]
    NonzeroCov <- which(Cov > 0)
    VCF        <- ReadVCF(VcfFileNames[i])
    if (sum(Cov) > 0){
      PeakSummary$FirstPos[i] <- min(NonzeroCov)
      PeakSummary$LastPos[i]  <- max(NonzeroCov)
      PeakSummary$MeanCovL1NonZero[i] <- mean(Cov[NonzeroCov])
      PeakSummary$MeanCovL1Range[i] <- mean(Cov[min(NonzeroCov):max(NonzeroCov)])
      PeakSummary$L1Width[i]  <- max(NonzeroCov) - min(NonzeroCov)
    }
    if (nrow(VCF) > 0){
      PeakSummary$wMismatch[i] <- sum(Cov[VCF$POS]) / sum(Cov)
    }
  }
  
  # Calculate correlation between different measures
  cor(PeakSummary[, c("LastPos","L1Width", "MaxCovL1", "MeanCovL1", "MeanCovL1NonZero", 
    "MeanCovL1Range", "MaxCovRef", "WidthRef", "wMismatch")],
    use = "pairwise.complete.obs")
  
  
  ######################################
  #                                    #
  #   Plots                            #
  #                                    #
  ######################################
  
  attach(PeakSummary) 
  plot(WidthRef, wMismatch, xlab = "Peak width [bp]", ylab = "Consensus mismatch")
  plot(WidthRef, L1Width, xlab = "Peak width [bp]", ylab = "L1 width [bp]")
  plot(WidthRef, MeanCovL1, xlab = "Peak width [bp]", ylab = "Mean L1 coverage")
  plot(MaxCovRef, L1Width, xlab = "Maximum peak coverage", ylab = "L1 width [bp]")
  plot(MaxCovRef, L1Width, xlab = "Maximum peak coverage", ylab = "L1 width [bp]", 
       xlim = c(0, 100))
  plot(MaxCovRef, L1Width, xlab = "Maximum peak coverage", ylab = "L1 width [bp]", 
       xlim = c(0, 50))
  plot(MaxCovRef, MeanCovL1, xlab = "Maximum peak coverage", ylab = "Mean L1 coverage", 
       xlim = c(0, 50))

  plot(MaxCovRef, MeanCovL1NonZero, xlab = "Maximum peak coverage", ylab = "Mean L1 coverage", 
       xlim = c(0, 50))
  plot(MaxCovRef, MeanCovL1Range, xlab = "Maximum peak coverage", ylab = "Mean L1 coverage", 
       xlim = c(0, 50))
  
  # pdf(file = CoverSummaryPlot)
  # plot(QuantileMat[2,], type = "n", ylim = c(0, max(QuantileMat)), 
  #      ylab = 'Coverage', xlab = "Genomic position")
  # polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
  #         col = "grey", border = NA)
  # lines(QuantileMat[2,], lwd = 1.2)
  # dev.off()
  
  # # Determine range index from file name 
  # pdf(file = CoverComparePlot)
  # plot(CoverMat[1,], type = "s", xlab = "Position on L1",
  #      ylab = "Coverage", ylim = c(0, 100))
  # Cols <- rainbow(nrow(CoverMat))
  # for (i in 1:nrow(CoverMat)){
  #   lines(CoverMat[i,], type = "s", col = Cols[i])
  #   
  # }
  # dev.off()
  
}


