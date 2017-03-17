##############################################
#
# General description:
#
#   The following function compares two types of genomic ranges with genomic
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
#    Freq1: frequency of GRanges3 names that are closest to GRanges1
#    Freq2: frequency of GRanges3 names that are closest to GRanges2
#    POverall: P-value for overall comparison 
#    PByName: vector of P-values for single comparisons
##############################################


Compare2RangesWith3rd <- function(GRanges1, GRanges2, GRanges3, Range3Names,
                                  Type12Names = c("Type1", "Type2")){
  
  # Test that GRanges3 and Range3Names have the same length
  if (length(GRanges3) != length(Range3Names)){
    stop("GRanges3 and Range3Names have to have the same length!\n")
  }
  
  # Find the nearest genomic range of GRanges3 for each, GRanges1 and GRanges1
  idxNearest1     <- nearest(GRanges1, GRanges3, ignore.strand = T)
  idxNearest2     <- nearest(GRanges2, GRanges3, ignore.strand = T)
  NearestNames1   <- as.character(Range3Names[idxNearest1])
  NearestNames2   <- as.character(Range3Names[idxNearest2])

  # Create matching vectors of closests 
  Range3NamesOrdered <- c(NearestNames1, NearestNames2)
  Type12NamesOrdered <- c(rep(Type12Names[1], length(idxNearest1)), 
                          rep(Type12Names[2], length(idxNearest2)))
  
  # Perform global chi-square tests
  ChisqTest <- chisq.test(Range3NamesOrdered, Type12NamesOrdered)
  
  # Test per segment for a difference in proportion
  Freq1    <- table(NearestNames1)
  Freq2    <- table(NearestNames2)
  Total1   <- sum(Freq1)
  Total2   <- sum(Freq2)
  
  # Determine which names are closest to both ranges 
  NamesInBoth <- intersect(NearestNames1, NearestNames2)
  NamesInBoth <- NamesInBoth[!is.na(NamesInBoth)]
  ORs <- Freq1[NamesInBoth]/Total1 / (Freq2[NamesInBoth]/Total2)
  x <- NamesInBoth[4]
  FisherTPVals   <- sapply(NamesInBoth, function(x){
    print(x)
    FTest <- fisher.test(Range3NamesOrdered == x, Type12NamesOrdered)
    FTest$p.value
  })
  
  # Create a output list
  list(Freq1 = Freq1, Freq2 = Freq2, ORs = ORs, ChisqTest = ChisqTest, 
       FisherTPVals = FisherTPVals, FisherTP_Bonf = FisherTPVals * length(FisherTPVals))
  
  
}
  