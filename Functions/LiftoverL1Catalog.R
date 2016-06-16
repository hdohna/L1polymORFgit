##############################################
#
# General description:
#
#   The following function lifts over L1 catalog coordinates from hg38
#   to hg19

# Input:
#
#     L1Catalog: table with L1 
#     ChainFilePath: path to chain file 

# Output (as list):
#   
#      NrMapped_hg19: per range the number of ranges it is mapped into
#      idxUniqueMapped: index of rows of original L1catalog that could be 
#         uniquely mapped
#      GRCatalogue_hg19: uniquely mapped genomic ranges on hg19


# Comments:
#    Requires bioconductor packages GenomicRanges and csaw (extractReads)

##############################################

LiftoverL1Catalog <- function(L1Catalog,
  ChainFilePath = "/home/hzudohna/L1polymORF/Data/hg38ToHg19.over.chain"){
  
  # Create genomic ranges of catalog on hg38
  GRCatalogue_hg38  <- GRanges(seqnames = L1Catalog$Chromosome,
                               ranges = IRanges(start = pmin(L1Catalog$start_HG38,
                                                             L1Catalog$end_HG38),
                                                end = pmax(L1Catalog$start_HG38,
                                                           L1Catalog$end_HG38)),
                               strand = L1Catalog$Strand)

  # Lift over genomic ranges to hg19
  GRCatalogue_hg19 <- liftOver(GRCatalogue_hg38, 
                               chain = import.chain(ChainFilePath))
  # Return output as list
  NrMapped_hg19    <- sapply(GRCatalogue_hg19, length)
  idxUniqueMapped  <- which(NrMapped_hg19 == 1) 
  list(NrMapped_hg19 = NrMapped_hg19,
      idxUniqueMapped = idxUniqueMapped ,
      GRCatalogue_hg19 = unlist(GRCatalogue_hg19[idxUniqueMapped]))
}

