##############################################
#
# General description:
#
#   The following function lifts genomic ranges over and only keeps uniquely
#   mapped ranges

# Input:
#
#     InputGRanges: genomic ranges to be lifted over
#     ChainFilePath: path to chain file 

# Output (as list):
#   
#      idxUniqueMapped: index of rows of original InputGRanges that could be 
#         uniquely mapped
#      LiftedRanges: GRanges that could be uniquely mapped to other genome built


# Comments:
#    Requires bioconductor packages GenomicRanges 

##############################################

UniqueLiftover <- function(InputGRanges,
  ChainFilePath = "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"){
  
  # Lift over genomic ranges 
  LiftedRangeList <- liftOver(InputGRanges, 
                               chain = import.chain(ChainFilePath))
  # Determine the unique mappings per read
  NrMapped    <- sapply(LiftedRangeList, length)
  idxUniqueMapped  <- which(NrMapped == 1) 
  
  # Get new genomic ranges 
  LiftedRanges <- unlist(LiftedRangeList[idxUniqueMapped])
  
  # Create list to return
  list(NrMapped = NrMapped,
      idxUniqueMapped = idxUniqueMapped,
      LiftedRangeList = LiftedRangeList,
      LiftedRanges = LiftedRanges)
}

