
##############################################
#
# General description:
#
#   The following function aggregates values that are associated with one set of
#   genomic ranges (GRanges2) within another set of genomci ranges (GRanges1)

# Input:
#
#     GRanges1: set of genomic ranges within values are aggregated
#     GRanges2: set of genomic ranges associated with the values
#     Values: vector of values (matching GRanges2) that are to be aggregated
#     within GRanges 1

# Output:
#   
#    AggData: 

##############################################


AggregateValsByGRanges <- function(GRanges1, GRanges2, Values,
                                   ValueName = "Mean",
                                 ValueFunction = mean){
  
  # Test that GRanges2 and Values have the same length
  if (length(GRanges2) != length(Values)){
    stop("Values have to match GRanges2!\n")
  }

  # Set up empty data frame
  AggData <- data.frame()

  # Aggregate values per GRanges1   
  Overlaps <- findOverlaps(GRanges1, GRanges2)
  if (length(Overlaps@from) > 0){
    AggData  <- aggregate(
      Values[Overlaps@from] ~ Overlaps@to, FUN = ValueFunction)
    colnames(AggData)  <- c("idxGR", ValueName)
  }
  AggData
}
  