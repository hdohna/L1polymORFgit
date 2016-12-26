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
  PosVector      <- Start:End
  PosNotNA       <- !is.na(ReadList$pos)
  ReadList$pos   <- ReadList$pos[PosNotNA]
  ReadList$cigar <- ReadList$cigar[PosNotNA]
  for (i in 1:length(ReadList$pos)){
    ReadStart <- ReadList$pos[i]
    ReadEnd   <- ReadStart + ReadLengthFromCigar(ReadList$cigar[i]) - 1
    CoverAdd <- (PosVector >= ReadStart) & (PosVector <= ReadEnd)
    CoverageVector <- CoverageVector + CoverAdd
  }
  CoverageVector
}



