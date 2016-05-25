# The script reads in L1 catalogue reference sequences, aligns consensus L1 and saves

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(ShortRead)
library(csaw)
library(seqinr)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue table and sequence 
L1CatalogTablePath <- "D:/L1polymORF/Data/L1Catalogue_Sat_May_07_15-15-31_2016.csv"
L1CatalogSeqPath   <- "D:/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fas"

# Create an output path
OutputPath <- gsub(".fas", "_L1Locations.RData", L1CatalogSeqPath)

# Minimum fragment size
MinFragSize <- 20

############################
#                          #
#    Read L1 catalogue     #
#                          #
############################

# Read in table with known L1 
L1Catalogue <- read.csv(L1CatalogTablePath, as.is = T)

# Create genomic 

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$L1Seq) 
blnL1Mapped2ref   <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

# Read in L1 catalog sequence file
L1withFlank <- read.fasta(L1CatalogSeqPath)
  
########################################
#                                      #
#    Get genomic ranges of all L1Hs    #
#      from repeatMasker               #
#                                      #
########################################

# Read repeat table
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]

# Subset to get only L1HS rows and save
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS",]

# Make some corrections to create a proper GRanges object with L1 Seqences
L1HSIRanges <- IRanges(start = RepeatTable$genoStart,
                       end = RepeatTable$genoEnd)
L1HSGRanges <- GRanges(seqnames = RepeatTable$genoName, ranges = L1HSIRanges,
                       strand = RepeatTable$strand)
L1HSGRanges <- L1HSGRanges[width(L1HSGRanges) < 5000]

#####################################
#                                   #
#    Analyze genomic ranges         #
#                                   #
#####################################

# Create genomic ranges from catalogue
L1CatalogueMapped <- L1Catalogue[!is.na(L1Catalogue$start_HG38) & L1Catalogue$Allele == 1,]
GRCatalogue <- GRanges(seqnames = L1CatalogueMapped$Chromosome,
                       ranges = IRanges(start = pmin(L1CatalogueMapped$start_HG38,
                                                     L1CatalogueMapped$end_HG38),
                                        end = pmax(L1CatalogueMapped$start_HG38,
                                                   L1CatalogueMapped$end_HG38)),
                       strand = L1CatalogueMapped$Strand)

# Get minimum distance between successive L1
GRprecede <- precede(GRCatalogue, GRCatalogue, ignore.strand = T)
min(abs(start(GRCatalogue) - start(GRCatalogue)[GRprecede]), na.rm = T)
cbind(as.vector(seqnames(GRCatalogue)), as.vector(seqnames(GRCatalogue))[GRprecede])

# Get all L1HS ranges that intersect with a 5kb flank around each L1
GRCatalogLarge <- resize(GRCatalogue, 16000, fix="center")

# Get indices of L1 genomic ranges that overlap with  
idxL1HSCatalogOverlaps <- findOverlaps(L1HSGRanges, GRCatalogLarge)

#####################################
#                                   #
#    Check for L1 in sequences      #
#                                   #
#####################################

# Get length of all L1s
L1withFlankLen <- sapply(L1withFlank, length)

# Get L1 consensus sequence
L1HSConsensus        <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1HSConsensusDNASt   <- DNAString(paste(L1HSConsensus[[1]], collapse = ""))
L1HSConsensusDNAStRC <- reverseComplement(L1HSConsensusDNASt)

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
names(AlignList) <- names(L1withFlank)

# Determine start and end of L1 in
L1StartEnd <- sapply(1:length(L1withFlank), function(x){
  c(Start = start(AlignList[[x]]@subject@range), 
    End = end(AlignList[[x]]@subject@range),
    Width = width(AlignList[[x]]@subject@range))
})
colnames(L1StartEnd) <- names(L1withFlank)

L1CatalogueNoL1 <- L1CatalogL1Mapped[L1StartEnd["Width", ] < 5500, ]
if (any(L1StartEnd["Width", ] < 5500)) {
  stop("Some catalogue sequences do not contain L1!\n")
} else {
  cat("L1 could be mapped to all", length(L1withFlank), "catalogue sequences\n")
  
}

# Create a list of local alignments of consensus sequence to catalogue
AlignListFlank <- lapply(1:length(L1withFlank), function(x){
  cat("Looking for L1 in seq", x, "of", length(L1withFlank), "\n")
  SubjectSeq <- DNAString(paste(L1withFlank[[x]], collapse = ""))
  SubjectSeqLeft  <- SubjectSeq[1:L1StartEnd[1,x]]
  SubjectSeqRight <- SubjectSeq[L1StartEnd[2,x]:length(SubjectSeq)]
  list(PwALeft = list(pos = pairwiseAlignment(L1HSConsensusDNASt, SubjectSeqLeft, 
                                              type = "local"),
                      neg = pairwiseAlignment(L1HSConsensusDNAStRC, SubjectSeqLeft, 
                                              type = "local")),
       PwARight = list(pos = pairwiseAlignment(L1HSConsensusDNASt, SubjectSeqRight, 
                                              type = "local"),
                      neg = pairwiseAlignment(L1HSConsensusDNAStRC, SubjectSeqRight, 
                                              type = "local")))                      
})
names(AlignListFlank) <- names(L1withFlank)

# Create ranges of L1 fragments
L1FragmentRanges <-lapply(1:length(AlignListFlank), function(i){
  IRs <- c(AlignListFlank[[i]]$PwALeft$pos@subject@range,
           AlignListFlank[[i]]$PwALeft$neg@subject@range,
           AlignListFlank[[i]]$PwARight$pos@subject@range,
           AlignListFlank[[i]]$PwARight$neg@subject@range)
  IRs <- IRs[width(IRs) > MinFragSize]
  if (length(IRs) > 0){
    GRanges(seqnames = names(L1withFlank)[i], ranges = IRs)
  }
})
NullElements <- sapply(L1FragmentRanges, is.null)
L1FragmentRanges <- GRangesList(L1FragmentRanges[!NullElements])
L1FragmentRanges <- unlist(L1FragmentRanges)

# Save AlignList and L1StartEnd
cat("Saving results in", OutputPath, "\n")
save(list = c("AlignList", "L1StartEnd", "L1FragmentRanges", "L1withFlank"), 
     file = OutputPath)
