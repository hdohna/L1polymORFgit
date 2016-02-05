# General description:
#
#    This function slides a query sequence along a subject sequence and finds
#    the offset that produces the best match between the two.

# Arguments:
#   
#    QuerySeq: vector with single letters
#    SubjectSeq: vector with single letters
#    HalfRange: integer values indicating half of the range of offset values to
#        be tested


# Output:
#   
#   AlignedQuery: vector of aligned query positions
#   AlignedSubject: vector of aligned subject positions 
#   OffSet: integer value giving best offset for query sequence
#   PropFit: proportion of query positions that fit with subject position
                

# Comment:

QuickAlign <- function(QuerySeq, SubjectSeq, HalfRange = 10){
  Range   <- (-HalfRange):HalfRange
  SeqPos  <- which(QuerySeq != "-")
  NrFit <- sapply(Range, function(x) {
    sum(SubjectSeq[SeqPos[(HalfRange + 1):(length(SeqPos) - HalfRange)]] == 
          QuerySeq[SeqPos[(HalfRange + 1 + x):(length(SeqPos) - HalfRange + x)]])})
  OffSet <- Range[which.max(NrFit)]
  Start1 <- max(1, 1 + OffSet)
  End1   <- min(length(SeqPos), length(SeqPos) + OffSet)
  Start2 <- max(1, 1 - OffSet)
  End2   <- min(length(SeqPos), length(SeqPos) - OffSet)
  list(
    AlignedQuery = QuerySeq[SeqPos[Start1:End1]],
    AlignedSubject = SubjectSeq[SeqPos[Start2:End2]],
    OffSet = OffSet,
    PropNotFit = sum(QuerySeq[SeqPos[Start1:End1]] != 
                    SubjectSeq[SeqPos[Start2:End2]]) / length(SeqPos)
  )
}



