# The following script creates a bed file for the L1 ranges for the reference
# genome hg19 with the non-reference L1 inserions

# Load packages
library(seqinr)
library(rtracklayer)

# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable   <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")
L1Length      <- RepeatTable$genoEnd - RepeatTable$genoStart
RepeatTable   <- RepeatTable[L1Length > 6000, ]
RepeatTableGR <- makeGRangesFromDataFrame(RepeatTable,
                         seqnames.field = "genoName",
                         start.field = "genoStart",
                         end.field = "genoEnd")


# Read in fasta file with the non-reference L1 ranges
L1withFlank <- read.fasta('D:/L1polymORF/Data/L1Catalog_NonRefWithFlank10000_Wed_Aug_10_17-32-20_2016.fas')

# Get L1 consensus sequence
L1HSConsensus        <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1HSConsensusDNASt   <- DNAString(paste(L1HSConsensus[[1]], collapse = ""))
L1HSConsensusDNAStRC <- reverseComplement(L1HSConsensusDNASt)

# Read in table with known L1 
L1Catalog <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", 
                        as.is = T)
AccMatch           <- match(names(L1withFlank), L1Catalog$Accession)
L1CatalogL1Matched <- L1Catalog[AccMatch, ]
blnL1Mapped        <- !is.na(L1Catalog$start_HG38)
blnL1Allele1       <- L1Catalog$Allele == 1
L1CatalogMapped    <- L1Catalog[blnL1Mapped & blnL1Allele1,]
L1Catalog_hg19     <- LiftoverL1Catalog(L1CatalogMapped, 
   ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")$L1CatalogWithHG19

# Create a genomic ranges object of the L1 catalog and match it with the repeat masker
L1CatalogGR <- makeGRangesFromDataFrame(L1Catalog_hg19,
                                          start.field = "start_HG19",
                                          end.field = "end_HG19",
                                          strand.field = "strand_L1toRef")
L1CatRptMsk_overlap <- findOverlaps(L1CatalogGR, RepeatTableGR, ignore.strand = T,
                                    minoverlap = 6000)

# Add an accession column to the repeat table
RepeatTable$Accession <- NA
RepeatTable$Accession[L1CatRptMsk_overlap@subjectHits] <- 
  L1Catalog_hg19$Accession[L1CatRptMsk_overlap@queryHits]

# Create a list of local alignments of consensus sequence to catalogue
AlignList <- lapply(1:length(L1withFlank), function(x){
  cat("Looking for L1 in seq", x, "of", length(L1withFlank), "\n")
  if (L1CatalogL1Matched$Strand[x] == "+"){
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

# Create a combined table and save it
L1Combined <- data.frame(chrom = c(as.character(RepeatTable$genoName), names(L1withFlank)),
                         start = c(RepeatTable$genoStart, StartsNonRef),
                         end   = c(RepeatTable$genoEnd, EndsNonRef),
                         name  = c(RepeatTable$Accession, names(L1withFlank)))
write.table(L1Combined, "D:/L1polymORF/Data/L1Ranges_hg19_withNonRefL1.txt",
            quote = F, col.names = T, row.names = F)


# L1NonRefGR   <- GRanges(seqnames = names(L1withFlank),
#                         ranges = IRanges(start = StartsNonRef,
#                                          end = EndsNonRef),
#                         strand = L1CatalogL1Matched$strand_L1toRef)
#   
#   
# # Create genomic ranges for L1s in reference and non-reference
# L1GR <- c(RepeatTableGR, L1NonRefGR)
# L1GR <- GRanges(seqnames = c(RepeatTable$genoName, names(L1withFlank)),
#   ranges = IRanges(start = c(RepeatTable$genoStart, StartsNonRef),
#                     end = c(RepeatTable$genoEnd, EndsNonRef)),
#   strand = c(as.character(RepeatTable$strand), L1CatalogL1Matched$strand_L1toRef))
L1GR <- makeGRangesFromDataFrame(L1Combined)
export.bed(L1GR, "D:/L1polymORF/Data/L1Ranges_hg19_withNonRefL1")
