# The script below reads the L1 

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Path to L1 catalogue file 
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Sat_May_07_15-15-31_2016.csv"

# Length of Flanking sequence to be used for alignment
FlankLength <- 300
TargetRange <- 100
ProductRange <- '400-500'

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
    SeqStart1  <- max(1, L1CatalogL1Mapped$start_Clone[i] - FlankLength)
    SeqEnd1    <- L1CatalogL1Mapped$start_Clone[i] + FlankLength
    SeqStart2  <- L1CatalogL1Mapped$end_Clone[i] - FlankLength
    SeqEnd2    <- min(length(SourceSeq), L1CatalogL1Mapped$end_Clone[i] + FlankLength)
    Seq1       <- toupper(SourceSeq[SeqStart1:SeqEnd1]) 
    Seq2       <- toupper(SourceSeq[SeqStart2:SeqEnd2]) 
    Seq1       <- paste(Seq1, collapse = "") 
    Seq2       <- paste(Seq2, collapse = "")  
    c(Seq1 = Seq1, Seq2 = Seq2)
  } else {
    Chrom      <- L1CatalogL1Mapped$Chromosome[i]
    SourceSeq  <- BSgenome.Hsapiens.UCSC.hg38[[Chrom]]
    SeqStart1  <- max(1, L1CatalogL1Mapped$start_HG38[i] - FlankLength)
    SeqEnd1    <- L1CatalogL1Mapped$start_HG38[i] + FlankLength
    SeqStart2  <- L1CatalogL1Mapped$end_HG38[i] - FlankLength
    SeqEnd2    <- min(length(SourceSeq), L1CatalogL1Mapped$end_HG38[i] + FlankLength)
    Seq1       <- SourceSeq[SeqStart1:SeqEnd1] 
    Seq2       <- SourceSeq[SeqStart2:SeqEnd2] 
    c(Seq1 = as.character(Seq1), Seq2 = as.character(Seq2))
  }
})
closebank()

#####################################
#                                   #
#    Write out Primer3 input        #
#                                   #
#####################################
  
# Lines for primer3 input
InputLineList <- lapply(1:length(L1withFlank), function(i){
  Seqs <- L1withFlank[[i]]
  c(paste('SEQUENCE_ID=', L1CatalogL1Mapped$Accession[i], "_S1", 
                         sep = ""),
  paste('SEQUENCE_TEMPLATE=', Seqs[1], sep = ""),
  paste('SEQUENCE_TARGET=', FlankLength, ",", TargetRange, sep = ""),
  paste('PRIMER_PRODUCT_SIZE_RANGE=', ProductRange, sep = ""),
  "PRIMER_TASK=generic", "=",
  paste('SEQUENCE_ID=', L1CatalogL1Mapped$Accession[i], "_S2", 
        sep = ""),
  paste('SEQUENCE_TEMPLATE=', Seqs[2], sep = ""),
  paste('SEQUENCE_TARGET=', FlankLength, ",", TargetRange, sep = ""),
  paste('PRIMER_PRODUCT_SIZE_RANGE=', ProductRange, sep = ""),
  "PRIMER_TASK=generic", "=")
})


# Write out L1 sequences with flank
PathSplit1 <- strsplit(L1CataloguePath, "L1Catalogue")[[1]]
PathSplit2 <- strsplit(PathSplit1[2], ".csv")[[1]][1]
OutputPath <- paste(PathSplit1[1], "P3input", PathSplit2, 
                    ".txt", sep = "")
con <- file(OutputPath, "wb")
writeLines(unlist(InputLineList), con = con)
close(con)


