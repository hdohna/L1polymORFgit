##############################################
#
# General description:
#
#   The following function reads a bam file and ranges of known L1 on the 
#   reference genome and returns ranges 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     ChmromLengths: named integer vector giving the length of each chromosome
#     L1RefRanges_All: ranges of all known L1 elements in reference genome
#     L1RefRanges_Funct: ranges of all known functional L1 elements in 
#       reference genome
#    MinMaxCover: integer giving minimum maximum coverage to be called a peak 
#    MinDist2L1: integer giving  minimum distance to L1 to be called a peak 

# Output:
#   
#    : ...

# Comments:
#    Requires bioconductor packages GenomicRanges and csaw (extractReads)

##############################################

CompareRanges <- function(
  BamFile           = "/home/hzudohna/BoData/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam",
  ChmromLengths,
  L1RefRanges_All,
  L1RefRanges_Funct,
  MinMaxCover = 40,
  MinGap = 6000,
  MinDist2L1  = 3*10^4){
  
  # Read coverage per chromosome
  CoverList <- lapply(1:length(ChmromLengths), function(i){
    Chrom <- names(ChmromLengths)[i]
    cat("Reading reads for chromosome", Chrom, "\n")
    ChromLength <- ChmromLengths[i]
    R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))
    Reads   <- extractReads(bam.file = BamFile , region = R1)
    ReadCov <- coverage(Reads)
  })
  
  #######################################
  #                                     #
  #    Determine 'islands' and          #
  #        overlap with L1              #
  #                                     #
  #######################################
  
  # Determine separate islands with continuous read coverage and turn islands 
  # into genomic ranges
  IslandList <- lapply(CoverList, function(x){
    Islands <- slice(x, lower = 1)
  })
  IslandGRanges <- lapply(1:length(IslandList), function(i){
    GRanges(seqnames = Chroms[i], 
            ranges = IslandList[[i]]@listData[[1]]@ranges,
            coverTotal = viewSums(IslandList[[i]])[[1]],
            coverMax   = viewMaxs(IslandList[[i]])[[1]])
  })
  IslandGRanges <- GRangesList(IslandGRanges)
  IslandGRanges <- unlist(IslandGRanges)
  
  # Merge ranges that are less than MinGap bp apart
  IslGRanges_reduced <- reduce(IslandGRanges, min.gapwidth = MinGap,
                               with.revmap = T)
  
  # Find overlaps between islands and L1 ranges
  blnOverlapIslands_All <- overlapsAny(IslGRanges_reduced, L1RefRanges_All)
  
  #######################################################
  #                                                     #
  #    Get ranges of suspected and known L1             #
  #                                                     #
  #######################################################
  
  # Determine maximum cover per island
  maxCoverOriginal   <- IslandGRanges@elementMetadata@listData$coverMax
  maxCover <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                     function(x) max(maxCoverOriginal[x]))
  
  # Find index of suspected L1 ranges (= islands with maximum coverage > 
  # MinMaxCover and that do not overlap with reference ranges)
  idxSuspectL1Ranges <- which(maxCover > MinMaxCover & (!blnOverlapIslands_All))
  SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
  
  # Remove ranges of suspected L1s that are too c
  DistToNearestL1    <- nearest(SuspectL1Ranges, L1GRanges)
  idxSuspectL1Ranges <- idxSuspectL1Ranges[DistToNearestL1 >= MinDist2L1]
  SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
  
  # Save results
  cat("*******  Saving results ...   *******\n")
  save(list = c("SuspectL1Ranges", "IslGRanges_reduced"), file = OutResults)
  
}

