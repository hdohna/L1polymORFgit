##############################################
#
# General description:
#
#   The following function process data of Tajima's D

# Input:
#
#     TajimaData: data.frame with Tajima's D data as produced by vcftools
#     L1GRanges:   GRanges of L1s
#     StepSize:    integer indicating the window step size
#     FlankSize:   integer indicating size of window around each L1
#     L1widthFull: integer indicating size of full-length L1

# Output:
#    blnFinished: boolean vector indicating for each job ID whether job has
#                    finished


# Comments:
#    

##############################################


ProcessTajimaDData <- function(TajimaData,
                               L1GRanges,
                               StepSize,
                               FlankSize,
                               L1widthFull = 6000){
  
  # Add a column with NA for intoerval
  TajimaData$TajimaD_modified <- TajimaData$TajimaD
  TajimaData$TajimaD_modified[TajimaData$N_SNPS == 0] <- NA
  
  # Intialize genomic ranges and matching dataframe with Tajima's D
  TajimaD_GR <- GRanges()
  TajimaD_reordered <- data.frame()
  
  # Loop over chromosomes, order positions per chromosome
  for (chr in unique(TajimaData$CHROM)){
    TDataSubset <- TajimaData[TajimaData$CHROM == chr,]
    startOrder  <- order(TDataSubset$BIN_START)
    TDataSubset <- TDataSubset[startOrder, ]
    New_GR <- GRanges(seqnames = chr, 
                      ranges = IRanges(start = TDataSubset$BIN_START[-nrow(TDataSubset)],
                                       end  = TDataSubset$BIN_START[-1]))  
    TDataSubset    <- TDataSubset[-1, ]
    
    # Retain only GRs and TData rows where New_GR have proper width
    blnProperWidth <- width(New_GR) == StepSize + 1
    New_GR         <- New_GR[blnProperWidth]
    TDataSubset    <- TDataSubset[blnProperWidth, ]
    
    # Append new data
    TajimaD_reordered <- rbind(TajimaD_reordered, TDataSubset)
    TajimaD_GR <- c(TajimaD_GR, New_GR)
  }
  TajimaD_GR <- as(TajimaD_GR, "GRanges")
  
  # Find ranges with only one L1
  TajDL1_OLcount <- countOverlaps(TajimaD_GR, L1GRanges)
  TajimaD_GR1    <- TajimaD_GR[TajDL1_OLcount == 1]
  
  # Get indices of L1 overlapping with TajimaD_GR1
  idxL1_OL <- which(overlapsAny(L1GRanges, TajimaD_GR1))
  
  # Get indices of L1s with non-overlapping neighborhoods
  idxUniqueL1Neighbor <- idxL1_OL[countOverlaps(L1Neighborhoods[idxL1_OL], 
                                                L1Neighborhoods[idxL1_OL]) == 1]
  
  # Find overlaps
  OL_L1TD <- findOverlaps(L1Neighborhoods[idxUniqueL1Neighbor], TajimaD_GR)
  
  # Create a list of Tajima's D values around each L1  
  TDList <- lapply(unique(OL_L1TD@from), function(x) {
    idxNeibhborL1 <- OL_L1TD@to[OL_L1TD@from == x]
    TDs <- TajimaD_reordered$TajimaD_modified[idxNeibhborL1]
    EdgeMean <- mean(TDs[c(1, length(TDs))], na.rm = T)
    if (is.na(EdgeMean)) {EdgeMean <- 0}
    TDs - EdgeMean
  })
  
  # Create a matrix of Tajima's D values
  NrInt <- FlankSize %/% StepSize
  MidPt <- NrInt %/% 2
  TDMat <- sapply(TDList, function(x){
    if (length(x) < NrInt){
      rep(NA, NrInt)
    } else {
      x[1:NrInt]
    }
  })
  
  # List of objects to output
  list(
    
    # Number of intervals and midpoint
    NrInt = NrInt,
    MidPt = MidPt,
    
    # matrix of Tajima's D values
    TDMat = TDMat,
    
    # Indices of fragment and full-length L1
    idxFull = which(width(L1GRanges[idxUniqueL1Neighbor][unique(OL_L1TD@from)]) >= L1widthFull),
    idxFragm = setdiff(1:length(unique(OL_L1TD@from)), idxFull),
    
    # Mean Tajima's D per fragment and full-length L1
    MeanTD_Full  = rowMeans(TDMat[,idxFull], na.rm = T),
    MeanTD_Fragm = rowMeans(TDMat[,idxFragm], na.rm = T),
    MedTD_Full   = apply(TDMat[,idxFull], 1, FUN = function(x) median(x, na.rm = T)),
    MedTD_Fragm  = apply(TDMat[,idxFragm], 1, FUN = function(x) median(x, na.rm = T))
    
  )
}