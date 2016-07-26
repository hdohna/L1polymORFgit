# The script below concatenates known L1 insertions and their surrounding into one reference 
# sequence

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(ShortRead)
library(csaw)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Sat_May_07_15-15-31_2016.csv"

# Length of Flanking sequence to be used for alignment
FlankLength <- 10000
FlankLengthShort <- 200

############################
#                          #
#    Read L1 catalogue     #
#                          #
############################

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$L1Seq) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[blnL1Mapped & blnAllele1,]

############################
#                          #
#    Obtain sequences      #
#                          #
############################

# Open connection to genbank
choosebank("genbank")

# Loop through sequences and obtain the flanking sequences
cat("Looping through sequences and obtain the flanking sequences\n\n")
L1Width <- L1CatalogL1Mapped$end_HG38 - L1CatalogL1Mapped$start_HG38
L1withFlank <- lapply(1:nrow(L1CatalogL1Mapped), function(i){
  
  cat("Analyzing entry", i, "of", nrow(L1CatalogL1Mapped), "\n")
  # Get sequence from clone if insertion is not present in reference
  if (is.na(L1CatalogL1Mapped$start_HG38[i]) | L1Width[i] < 6000){
    AccNr     <- L1CatalogL1Mapped$Accession[i]
    x         <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
    SourceSeq <- getSequence(x$req)[[1]]
    SeqStart  <- max(1, L1CatalogL1Mapped$start_Clone[i] - FlankLength)
    SeqEnd    <- min(length(SourceSeq), L1CatalogL1Mapped$end_Clone[i] + FlankLength)
    Seq       <- SourceSeq[SeqStart:SeqEnd] 
    Seq       <- toupper(Seq)
  } else {
    Chrom     <- L1CatalogL1Mapped$Chromosome[i]
    SourceSeq <- BSgenome.Hsapiens.UCSC.hg38[[Chrom]]
    SeqStart  <- max(1, L1CatalogL1Mapped$start_HG38[i] - FlankLength)
    SeqEnd    <- min(length(SourceSeq), L1CatalogL1Mapped$end_HG38[i] + FlankLength)
    Seq       <- SourceSeq[SeqStart:SeqEnd] 
    Seq       <- s2c(as.character(Seq))
  }
})
closebank()

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

# Determine start and end of L1 in
L1StartEnd <- sapply(1:length(L1withFlank), function(x){
  c(Start = start(AlignList[[x]]@subject@range), 
    End = end(AlignList[[x]]@subject@range),
    Width = width(AlignList[[x]]@subject@range))
})

L1CatalogueNoL1 <- L1CatalogL1Mapped[L1StartEnd["Width", ] < 5500, ]
if (any(L1StartEnd["Width", ] < 5500)) {
  stop("Some catalogue sequences do not contain L1!\n")
} else {
  cat("L1 could be mapped to all", length(L1withFlank), "catalogue sequences\n")
  
}

# # Align consensus to all L1 sequences
# AlignConsens2L1 <- lapply(1:nrow(L1CatalogL1Mapped), function(x){
#   cat("Aligning L1 consensus to seq", x, "of", nrow(L1CatalogL1Mapped), "\n")
#   pairwiseAlignment(L1HSConsensusDNASt, L1CatalogL1Mapped$L1Seq[x], 
#                     type = "global")
# })
# 
# # Count indels and mismatches 
# MismatchCount <- sapply(AlignConsens2L1, function(x){
#   c(MisMatch     = length(x@pattern@mismatch[[1]]),
#     IndelPattern = sum(width(x@pattern@indel)[[1]]),
#     IndelSubject = sum(width(x@subject@indel)[[1]]))
# })
# hist(colSums(MismatchCount))
# 
# # Check whether 
# blnL1present <- L1StartEnd["Width", ] > 5500
# boxplot(colSums(MismatchCount) ~ blnL1present)
# boxplot(MismatchCount["MisMatch",] ~ blnL1present)
# boxplot(MismatchCount["IndelPattern",] ~ blnL1present)
# boxplot(MismatchCount["IndelSubject",] ~ blnL1present)
# 
# # Create a list of local alignments of L1 sequence to catalogue
# AlignListSelf <- lapply(1:length(L1withFlank), function(x){
#   cat("Looking for L1 in seq", x, "of", length(L1withFlank), "\n")
#   if (L1CatalogL1Mapped$Strand[x] == "+"){
#     PatternSeq <- L1CatalogL1Mapped$L1Seq[x]
#   } else {
#     PatternSeq <- DNAString(L1CatalogL1Mapped$L1Seq[x])
#     PatternSeq <- reverseComplement(PatternSeq)
#   }
#   SubjectSeq <- DNAString(paste(L1withFlank[[x]], collapse = ""))
#   pairwiseAlignment(PatternSeq, SubjectSeq, type = "global")
# })
# 
# sapply(1:length(L1withFlank), function(x){
#   c(Start = start(AlignListSelf[[x]]@subject@range), 
#     End = end(AlignListSelf[[x]]@subject@range),
#     Width = width(AlignListSelf[[x]]@subject@range))
# })

#####################################
#                                   #
#    Write out L1 with flanks       #
#                                   #
#####################################
  
# Write out L1 sequences with flank
PathSplit1 <- strsplit(L1CataloguePath, "L1Catalogue")[[1]]
PathSplit2 <- strsplit(PathSplit1[2], ".csv")[[1]][1]
OutputPath <- paste(PathSplit1[1], "L1CatalogueWithFlank", FlankLength,
                    PathSplit2, 
                    ".fas", sep = "")
write.fasta(L1withFlank, L1CatalogL1Mapped$Accession, 
            file.out = OutputPath)


