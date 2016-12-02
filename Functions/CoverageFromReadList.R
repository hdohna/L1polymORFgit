# General description:
#
#    This function calculates the coverage from

# Arguments:
#   
#    ReadList: list of reads as produced by Rsamtools::scanBam
#    Start: integer giving starting position
#    End: integer giving end position

# Output:
#   
#    vector with coverage

# Comment: requires function ReadLengthFromCigar

CoverageFromReadList <- function(ReadList, Start = 1, End){
  CoverageVector <- rep(0, End - Start + 1)
  PosVector  <- Start:End
  ReadList$pos   <- ReadList$pos[!is.na(ReadList$pos)]
  ReadList$cigar <- ReadList$cigar[!is.na(ReadList$pos)]
  ReadStarts <- ReadList$pos
  ReadEnds   <- ReadStarts + sapply(ReadList$cigar, ReadLengthFromCigar) - 1
  for (i in 1:length(ReadList$pos)){
    ReadStart <- ReadList$pos[i]
    ReadEnd   <- ReadStart + ReadLengthFromCigar(ReadList$cigar) - 1
    CoverAdd <- (PosVector >= ReadStart) & (PosVector <= ReadEnd)
    CoverageVector <- CoverageVector + CoverAdd
  }
  CoverageVector
}



