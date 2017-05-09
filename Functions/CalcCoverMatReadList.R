##############################################
#
# General description:
#
#   The following function  reads in reads aligned to L1, creates a list of reads, 
#   calculates a coverage and quantile matrix
#   

# Input:
#
#     OutputFilePath: (character) path to file that will be written out
#     OutFolderName_NonRef: (character) path to folder that contains bam file
#         per peak
#     L1RangesPath: (character) path to file that contains results of comparing
#         peaks with L1 ranges (created by function ComparePeaksWithRefL1)

# Output:
#   
#     ...

##############################################

CalcCoverMatReadList <- function(
                         OutputFilePath,
                         OutFolderName_NonRef,
                         L1RangesPath,
                         GenomeBamPath){
  
  cat("***********************************************************\n")
  cat("**                                                       **\n")
  cat("**    Running function CalcCoverMatReadList ...          **\n")
  cat("**                                                       **\n")
  cat("***********************************************************\n")
  
  # Load ranges
  load(L1RangesPath)
  
  # get names of newly created bam files
  FileNames <- list.files(OutFolderName_NonRef, pattern = ".bam",
                          full.names = T)
  FileNames <- FileNames[-grep(".bam.", FileNames)]
  
  # Loop through file names and read in bam files of reads mapped to L1
  cat("*******  Scanning bam files per L1   *************\n")
  ScannedL1Ranges <- lapply(1:length(FileNames), function(i) {
    cat("Scanning bam file", i, "of", length(FileNames), "\n")
    scanBam(FileNames[i])
  })
  
  # Count the number of reads mapped
  NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
    sum(!is.na(x[[1]]$pos))
  })
  idxFilesWithReads <- which(NrMapped2L1 > 0)
  FilesWithReads    <- FileNames[idxFilesWithReads]
  
  # Extract the names of files with reads
  PeakNames <- sapply(FilesWithReads, function(x) {
    FPathSplit <- strsplit(x, "/")[[1]]
    FName <- FPathSplit[length(FPathSplit)]
    substr(FName, 1, nchar(FName) - 4)
  })
  
  # Get read list per peak
  ReadListPerPeak <- lapply(idxFilesWithReads, function(x) {
    ScannedL1Ranges[[x]][[1]]
  })
  names(ReadListPerPeak) <- PeakNames
  
  # Create a list that matches reads mapped to L1 and reads mapped
  # to the genome
  cat("*******  Calculating MatchedReadList   *************\n")
  MatchedReadList <- lapply(1:length(ReadListPerPeak), function(x) {
    
    # Extract range index from file name
    FName    <- names(ReadListPerPeak)[x]
    idxRange <- as.numeric(strsplit(FName, "_")[[1]][2])
    PeakGR    <- IslGRanges_reduced[idxRange]
    
    # Get reads mapped to L1, retain only primary reads and get the number of bp
    # clipped on the left
    RL          <- ReadListPerPeak[[x]]
    primMap     <- RL$flag <= 2047 & !(is.na(RL$pos))
    RL          <- lapply(RL, function(y) y[primMap])
    LRClippedL1 <- sapply(RL$cigar, NrClippedFromCigar, USE.NAMES = F)
    L1Length    <- width(RL$seq) - LRClippedL1[1,]
    
    # Get reads mapped to the genome on current locus
    ScanParam <- ScanBamParam(what = scanBamWhat(), which = PeakGR)
    GenomeRL  <- scanBam(GenomeBamPath,  param = ScanParam)
    # primMap     <- GenomeRL[[1]]$flag <= 2047 & !(is.na(GenomeRL[[1]]$pos))
    # GenomeRL[[1]]    <- lapply(GenomeRL[[1]], function(y) y[primMap])
    LRClippedG <- sapply(RL$cigar, NrClippedFromCigar)
    ReadMatch <- match(RL$qname, GenomeRL[[1]]$qname)
    if (any(is.na(ReadMatch))){
      browser()
    }
    
    # Generate a list of reads mapped to genome, matching the list of reads 
    # mapped to L1
    GenomeMatchRL <- lapply(GenomeRL[[1]], function(y) y[ReadMatch])
    list(L1ReadList = RL, GenomeReadList = GenomeMatchRL)
  })
  names(MatchedReadList) <- names(ReadListPerPeak)
  
  
  cat("*******  Calculating coverage matrix   *************\n")
  CoverMat <- t(sapply(1:length(ReadListPerPeak), function(i) {
    cat("Analyzing row", i, "of", length(ReadListPerPeak), "\n")
    x <- ReadListPerPeak[[i]]
    primMap <- x$flag <= 2047
    RL <- lapply(x, function(y) y[primMap])
    CoverageFromReadList(RL, End = 6064)
  }))
  dim(CoverMat)
  rownames(CoverMat) <- names(ReadListPerPeak)
  
  # Get means and 95% quantiles 
  cat("*******  Calculating quantile matrix   *************\n")
  QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
  idxFw <- 1:ncol(CoverMat)
  idxRv <- ncol(CoverMat):1
  plot(QuantileMat[2,], type = "n", ylim = c(0, max(QuantileMat)), 
       ylab = 'Coverage', xlab = "Genomic position")
  polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
          col = "grey", border = NA)
  lines(QuantileMat[2,], lwd = 1.2)
  
  # Save objects to file
  cat("*****  Saving data to", OutputFilePath,  "*********\n")
  save(list = c("ScannedL1Ranges", "idxFilesWithReads", "ReadListPerPeak", 
                "MatchedReadList","CoverMat", "QuantileMat", "FilesWithReads",
                "IslGRanges_reduced", "idxSuspectL1Ranges"), 
       file = OutputFilePath)
}