# General description:
#
#    This function converts starting positions of reads into genomic
#    ranges.
#    

# Arguments:
#   
#    BamFile: character string for path to bamfile

# Output:
#   
#    GRanges object for starting positions of reads

# Comment: requires function scanBam from Rsamtools

Reads2GRanges <- function(BamFile){
  ReadList <- scanBam(BamFile)[[1]]
  GRanges(seqnames = ReadList$rname, IRanges(start = ReadList$pos,
                                             end = ReadList$pos))
}