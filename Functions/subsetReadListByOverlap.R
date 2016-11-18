# General description:
#
#    This function subsets a list of reads (created by scanBam) to retain
#    only entries of reads that overlap with a particular genomic range
#
# Arguments:
#   
#    ReadList: list of reads (created by scanBam)
#    GR:  genomic range

# Output:
#   
#    ReadListSubset: subsetted read list

# Comment:

subsetReadListByOverlap <- function(ReadList, GR){
  
  ReadListSubset <- lapply(ReadList, function(LocalReadList){
    ReadGR <- GRanges(seqnames = LocalReadList$rname[1],
                      ranges = IRanges(LocalReadList$pos, 
                                       width = width(LocalReadList$seq)))
    blnOverlap <- overlapsAny(ReadGR, GR)
    lapply(LocalReadList, function(x) x[blnOverlap])
  })
}



