##############################################
#
# General description:
#
#   The following function aggregates values for two types of genomic ranges with genomic
#   ranges of multiple types and tests whether any type is enriched in one of
#   the two sets  

# Input:
#
#     GRanges1: set of genomic ranges of type 1
#     GRanges2: set of genomic ranges of type 2
#     GRanges3: set of genomic ranges of type 3
#     Range3Names: names of genomic range of type 3
#     Type12Names: names of GRanges1 and GRanges2, respectively

# Output:
#   
#    AggData: Dataframe that first contains aggregated value for GRanges1 and
#       then for ranges 2. Two columns: Values and Types

##############################################


AggregateValsBy2GRangesSet <- function(GRanges1, GRanges2, GRanges3, Values,
   Type12Names = c("full", "fragm"), ValueName = "Value", TypeName = "Type"){
  
  # Test that GRanges3 and Range3Names have the same length
  if (length(GRanges3) != length(Values)){
    stop("GRanges3 and Range3Names have to have the same length!\n")
  }
  
  # Set up empty data frames
  AggData1 <- data.frame()
  AggData2 <- data.frame()
  
  # Aggregate values per GRanges1   
  Overlaps1  <- findOverlaps(GRanges3, GRanges1)
  if (length(overlapL1Fragm@from) > 0){
    AggData1  <- aggregate(
      Values[Overlaps1@from] ~ Overlaps1@to, FUN = mean)
    colnames(AggData1)  <- c("idxGR", ValueName)
    AggData1$start      <- start(GRanges1)[AggData1$idxGR]
    AggData1$end        <- end(GRanges1)[AggData1$idxGR]
    AggData1$chromosome <- as.vector(seqnames(GRanges1))[AggData1$idxGR]
    AggData1$Type       <- Type12Names[1]
    colnames(AggData1)[colnames(AggData1) == "Type"] <- TypeName
  }
  
  # Aggregate histone values per full-length L1   
  Overlaps2  <- findOverlaps(GRanges3, GRanges2)
  if (length(Overlaps2@from) > 0){
    AggData2  <- aggregate(
    Values[Overlaps2@from] ~ Overlaps2@to, FUN = mean)
    colnames(AggData2)  <- c("idxGR", ValueName)
    AggData2$start      <- start(GRanges2)[AggData2$idxGR]
    AggData2$end        <- end(GRanges2)[AggData2$idxGR]
    AggData2$chromosome <- as.vector(seqnames(GRanges2))[AggData2$idxGR]
    AggData2$Type       <- Type12Names[2]
    colnames(AggData2)[colnames(AggData2) == "Type"] <- TypeName
  }
  # Combine both 
  rbind(AggData1, AggData2)
  
}
  