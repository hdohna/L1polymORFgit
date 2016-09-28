# The following script creates a bed file for the L1 ranges for the reference
# genome hg19 with the non-reference L1 inserions

# Load packages
library(seqinr)
library(rtracklayer)

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.csv("D:/L1polymORF/Data/repeatsHg38_L1HS.csv")
L1Length <- RepeatTable$genoEnd - RepeatTable$genoStart
RepeatTable <- RepeatTable[L1Length > 6000, ]

# Read in fasta file with the non-reference L1 ranges
L1withFlank <- read.fasta('D:/L1polymORF/Data/L1Catalog_NonRefWithFlank10000_Wed_Aug_10_17-32-20_2016.fas')

# Get L1 consensus sequence
L1HSConsensus        <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1HSConsensusDNASt   <- DNAString(paste(L1HSConsensus[[1]], collapse = ""))
L1HSConsensusDNAStRC <- reverseComplement(L1HSConsensusDNASt)

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", 
                        as.is = T)
AccMatch <- match(names(L1withFlank), L1Catalog$Accession)
L1CatalogL1Mapped <- L1Catalog[AccMatch, ]

# Create a list of local alignments of consensus sequence to catalogue
AlignList <- lapply(1:length(L1withFlank), function(x){
  cat("Looking for L1 in seq", x, "of", length(L1withFlank), "\n")
  if (L1CatalogL1Mapped$Strand[x] == "+"){
    PatternSeq <- L1HSConsensusDNASt
  } else {
    PatternSeq <- L1HSConsensusDNAStRC
  }
  SubjectSeq <- DNAString(paste(L1withFlank[[x]], collapse = ""))
  pairwiseAlignment(PatternSeq, SubjectSeq, type = "local")
})

# Get Starts and ends
StartsNonRef <- sapply(AlignList, function (x) start(x@subject@range))
EndsNonRef   <- sapply(AlignList, function (x) end(x@subject@range))

# Create genomic ranges for L1s in reference and non-reference
L1GR <- GRanges(seqnames = c(RepeatTable$genoName, names(L1withFlank)),
  ranges = IRanges(start = c(RepeatTable$genoStart, StartsNonRef),
                    end = c(RepeatTable$genoEnd, EndsNonRef)),
  strand = c(as.character(RepeatTable$strand), L1CatalogL1Mapped$strand_L1toRef))
export.bed(L1GR, "D:/L1polymORF/Data/L1Ranges_hg19_withNonRefL1")
