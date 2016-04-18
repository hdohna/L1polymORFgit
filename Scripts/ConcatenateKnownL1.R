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

# Length of Flanking sequence to be used for alignment
FlankLength <- 5000
blnCreateAlignment <- F

############################
#                          #
#    Read L1 catalogue     #
#                          #
############################

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalogue_Fri_Apr_01_18-28-08_2016.csv", 
                        as.is = T)

# Retain only entries with mapped L1 insertions
idxL1Mapped       <- !is.na(L1Catalogue$L1Seq) 
L1CatalogL1Mapped <- L1Catalogue[idxL1Mapped,]

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
    SeqStart  <- max(1, L1CatalogL1Mapped$start_Clone - FlankLength)
    SeqEnd    <- min(length(SourceSeq), L1CatalogL1Mapped$end_Clone + FlankLength)
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

# Write out L1 sequences with flank
write.fasta(L1withFlank, L1CatalogL1Mapped$Accession, 
            file.out = "D:/L1polymORF/Data/L1CatalogueWithFlank.fas")


############################
#                          #
#    Align sequences       #
#                          #
############################

if (blnCreateAlignment){
  
  # Write sequences as fasta file and align
  L1Seqs <- lapply(L1CatalogL1Mapped$L1Seq, function(x){
    s2c(x)
  })
  write.fasta(L1Seqs, L1Catalogue$Accession, 
              file.out  = "D:/L1polymORF/Data/L1CatalogueSeqsUnaligned.fas")
  run_MUSCLE(InputPath  = "D:/L1polymORF/Data/L1CatalogueSeqsUnaligned.fas", 
             OutputPath = "D:/L1polymORF/Data/L1CatalogueSeqsAligned.fas")
  
}

