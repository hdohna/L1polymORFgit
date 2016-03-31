##############################################
#
# General description:
#
#   The following script reads in various data sources ona repeat table downloaded from the genome
#   browser repeatMasker track (http://genome.ucsc.edu/cgi-bin/hgTables)
#   and subsets to get all L1HS ranges

# Input:
#
#    D:/L1polymORF/Data/repeatsHg38: table with all repeats
#   

# Output:
#   
#    L1HS_repeat_table.csv: csv file with all L1HS repeats
#    L1Sequences_reptab.fas: fasta file with all L1 sequences

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Length of sequences flanking L1 insertions [bp] used to locate insertions
FlankSize <- 100

# Boolean indicators for whether to perform particular processes
blnBuildBrouha2003 <- T
blnBuildBeck2010   <- F
blnFindVarSites    <- F
blnAlignsequences  <- F

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Load consensus sequence
L1Consens <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1Consens <- paste(L1Consens[[1]], collapse = "")
L1Consens <- DNAString(L1Consens)

# Read in repeatMasker tables 
L1repMask_Hg19   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv")
L1repMask_Hg38   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table.csv")

#######################################
#                                     #
#     Construct GRanges objects       #
#                                     #
#######################################

# Function to construct genomic ranges from data table
GRangesFromTable <- function(DataTable, SeqNCol, StartCol, EndCol, StrandCol = NULL){
  NewStart <- pmin(DataTable[,StartCol], DataTable[,EndCol])
  NewEnd   <- pmax(DataTable[,StartCol], DataTable[,EndCol])
  if (!is.null(StrandCol)){
    Strand   <- DataTable[,StrandCol]
  } else {
    Strand   <- NULL
  }
  idxNoNA  <- (!is.na(NewStart)) & (!is.na(NewEnd))
  IR <- IRanges(start = NewStart[idxNoNA], end = NewEnd[idxNoNA])
  GRanges(seqnames = DataTable[idxNoNA, SeqNCol], ranges = IR, strand = Strand)
}

# Genomic ranges for each dataset
GRanges_L1repMask_Hg19 <- GRangesFromTable(L1repMask_Hg19, "genoName", "genoStart", 
                                           "genoEnd", "strand")
GRanges_L1repMask_Hg38 <- GRangesFromTable(L1repMask_Hg38, "genoName", "genoStart", 
                                      "genoEnd", "strand")

# Function to remove genomic ranges below a lower limit and above an upper limit
RemoveRanges <- function(Ranges, LowerW = 0, UpperW = Inf){
  Ranges[width(Ranges) >= LowerW & width(Ranges) <= UpperW]
}

#######################################
#                                     #
#     Get flanking sequence of        #
#         full-length L1              #
#                                     #
#######################################

# Remove strange regions from GRanges_L1repMask_Hg38
GRanges_L1repMask_Hg38 <- RemoveRanges(GRanges_L1repMask_Hg38, LowerW = 6000, 
                                       UpperW = 12000)

# Get the 500 nuc upstream and downstream of the L1
Flank5PSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                     names  = seqnames(GRanges_L1repMask_Hg38), 
                     start  = start(GRanges_L1repMask_Hg38)  - FlankSize, 
                     end    = start(GRanges_L1repMask_Hg38), 
                     strand = strand(GRanges_L1repMask_Hg38)) 
Flank3PSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                     names  = seqnames(GRanges_L1repMask_Hg38), 
                     start  = end(GRanges_L1repMask_Hg38), 
                     end    = end(GRanges_L1repMask_Hg38) + FlankSize, 
                     strand = strand(GRanges_L1repMask_Hg38)) 


L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges_L1repMask_Hg38)

#################################################
#                                               #
#   Process L1 table form Brouha et al. 2003    #
#                                               #
#################################################

if (blnBuildBrouha2003){
  
  # Read in table of hot L1 (obtained from Brouha et al.2003 PNAS)
  cat("Read table with hot L1s \n")
  L1Brouha2003_raw <- read.delim('D:/L1polymORF/Data/L1HotTable_raw.txt', 
                                 sep = " ", as.is = T)

  # Find duplicated column names 
  CNames        <- colnames(L1Brouha2003_raw)
  ColNameSuffix <- substr(CNames, nchar(CNames) - 1, nchar(CNames))
  blnDuplCols   <- ColNameSuffix == ".1"
  
  # Separate table according to duplicated column names and append duplucated part
  L1Brouha2003Dupl <- L1Brouha2003_raw[,blnDuplCols]
  colnames(L1Brouha2003Dupl) <- CNames[!blnDuplCols]
  L1Brouha2003Table    <- rbind(L1Brouha2003_raw[,!blnDuplCols],L1Brouha2003Dupl)
  
  # Add chromosome column
  LChar          <- as.character(L1Brouha2003Table$Locus)
  LocusSpl1      <- sapply(LChar, function(x) strsplit(x, "[p-q]")[[1]][1])
  L1Brouha2003Table$Chr <- paste("chr", LocusSpl1, sep = "")
  
  # Rename columns
  colnames(L1Brouha2003Table)[colnames(L1Brouha2003Table) == "Accession_no."] <- "Accession"
  
  # Add more columns
  L1Brouha2003Table$L1Seq            <- NA
  L1Brouha2003Table$L1SeqFlank5p     <- NA
  L1Brouha2003Table$L1SeqFlank3p     <- NA
  L1Brouha2003Table$L1SeqSourceType  <- "BAC"
  L1Brouha2003Table$start_HG38       <- NA
  L1Brouha2003Table$end_HG38         <- NA
  L1Brouha2003Table$Strand           <- NA
  L1Brouha2003Table$Reference        <- "Brouha2003"
  
  # Download sequences 
  cat("Get sequences for hot L1 loci \n")
  load("D:/L1polymORF/Data/Brouha2003Sequences.RData")

  # Match each L1HS to a locus 
  cat("Match hot L1 loci to L1s from repeatMasker\n")
  L1SeqHot <- toupper(L1SeqHot)
  
  # Get strand and chromosome from repeatMasker table
  Strand_L1repMask_Hg38 <- as.vector(strand(GRanges_L1repMask_Hg38))
  chr_L1repMask_Hg38    <- as.vector(seqnames(GRanges_L1repMask_Hg38))
  
  # Find sequences that are matched by flanks
  L1SeqHotDNAStSet   <- DNAStringSet(L1SeqHot)
  L1SeqHotDNAStSetRV <- reverseComplement(L1SeqHotDNAStSet)
  
  # Loop through flanking sequences and find parts of the BAC clones that match
  # flanking sequences
  for (i in 1:length(Flank5PSeq)){
    print(i)
    if (Strand_L1repMask_Hg38[i] == "-"){
      DNAStSet <- L1SeqHotDNAStSetRV
      LeftP    <- Flank3PSeq[[i]]
      RightP   <- Flank5PSeq[[i]]
    } else {
      DNAStSet <- L1SeqHotDNAStSet
      LeftP    <- Flank5PSeq[[i]]
      RightP   <- Flank3PSeq[[i]]
    }
    idxChr   <- which(L1Brouha2003Table$Chr == chr_L1repMask_Hg38[i])
    DNAStSet <- DNAStSet[idxChr]
    FlankMatches <- lapply(DNAStSet, function(x){
      matchLRPatterns(LeftP, RightP, max.gaplength = 6500, 
                      subject = x, max.Lmismatch = 5, max.Rmismatch = 5,
                      with.Lindels = T, with.Rindels = T)
    })
    idxMatch <- which(sapply(FlankMatches, length) > 0)
    if (length(idxMatch) > 0) {
      mViews   <- FlankMatches[idxMatch][[1]]
      idxAbove6000 <- width(mViews) > 6000 & width(mViews) < 7000
      if (sum(idxAbove6000) == 1){
        S <- start(mViews[idxAbove6000])
        E <- end(mViews[idxAbove6000])
        Seq <- mViews[idxAbove6000]@subject[S:E] 
        Seq5P <- mViews[idxAbove6000]@subject[(S - FlankSize):(S - 1)] 
        Seq3P <- mViews[idxAbove6000]@subject[(E + 1):(E + FlankSize)] 
        L1Brouha2003Table$L1SeqFlank5p[idxChr[idxMatch]]  <- as.character(Seq5P)
        L1Brouha2003Table$L1SeqFlank3p[idxChr[idxMatch]]  <- as.character(Seq3P)
        L1Brouha2003Table$L1Seq[idxChr[idxMatch]] <- as.character(Seq)
        L1Brouha2003Table$start_HG38[idxChr[idxMatch]] <- start(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$end_HG38[idxChr[idxMatch]]   <- end(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$Strand[idxChr[idxMatch]] <- Strand_L1repMask_Hg38[i]
      } else {
        cat("More than one set of flanking sequences enclose a stretch above 6000\n")
      }
      
    }
  }
  
  # Break L1Brouha2003Table inot entries with and without L1 insertion on
  # the reference genome
  L1Brouha2003Ref    <- L1Brouha2003Table[!is.na(L1Brouha2003Table$L1Seq),]
  L1Brouha2003NonRef <- L1Brouha2003Table[is.na(L1Brouha2003Table$L1Seq),]
  
  # Get accession numbers and sequences that have no matching flanking sequence
  # in reference genome
  AccNonRef   <- L1Brouha2003NonRef$Accession
  L1SeqNonRef <- L1SeqHot[match(AccNonRef, names(L1SeqHot))] 
  
  # Loop through sequences that have no coordinates in the reference genome
  # and get the sequences through pairwise alignment
  L1DFBrouha2003  <- getNonRefL1(L1Consens, AccNrs = AccNonRef, 
                                 Seqs = DNAStringSet(L1SeqNonRef),
                               MinMatchWidth = 5500, FlankSize = FlankSize)
  
  # Replace each colum
  for (Col in colnames(L1DFBrouha2003)){
    L1Brouha2003NonRef[,Col] <- L1DFBrouha2003[,Col]
  }
  
  # Merge both table
  L1Brouha2003Merged <- merge(L1Brouha2003Ref, L1Brouha2003NonRef, all = T)
  
  # Write out table with sequence data
  write.csv(L1Brouha2003Merged, "D:/L1polymORF/Data/L1Brouha2003.csv")
  
} else {
    L1Brouha2003Merged <- read.csv("D:/L1polymORF/Data/L1Brouha2003.csv", 
                                   as.is = T)
}

#################################################
#                                               #
#   Process L1 table from Beck et al. 2003      #
#                                               #
#################################################

if (blnBuildBeck2010) {
  Beck2010Table <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withEmptySite.csv")
  L1DFBeck2010  <- getNonRefL1(L1Consens, AccNrs = Beck2010Table$Accession, 
                               MinMatchWidth = 5500, FlankSize = FlankSize)
  Beck2010TableWithL1 <- cbind(Beck2010Table, L1DFBeck2010)  
  write.csv(Beck2010TableWithL1, "D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv")
  
} else {
  Beck2010TableWithL1 <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv",
                                  as.is = T)
}

# # Create a dataframe that keeps track of various measures that help locate L1 
# # insertion on reference genome
# InsertDF <- data.frame(FlankStart        = rep(NA, nrow(Beck2010TableWithL1)),
#                        InsertionStartRel = rep(NA, nrow(Beck2010TableWithL1)),
#                        InsertionStartAbs = rep(NA, nrow(Beck2010TableWithL1)),
#                        NrNuc5p           = rep(NA, nrow(Beck2010TableWithL1)),
#                        NrNuc3p           = rep(NA, nrow(Beck2010TableWithL1)))
# 
# InsertDF <- data.frame(FlankStart        = rep(NA, nrow(L1Brouha2003Merged)),
#                        InsertionStartRel = rep(NA, nrow(L1Brouha2003Merged)),
#                        InsertionStartAbs = rep(NA, nrow(L1Brouha2003Merged)),
#                        NrNuc5p           = rep(NA, nrow(L1Brouha2003Merged)),
#                        NrNuc3p           = rep(NA, nrow(L1Brouha2003Merged)))
# 
# # Get the total length of the flanking sequence
# FlankTotal <- nchar(L1Brouha2003Merged$L1SeqFlank5p2x) + 
#   nchar(L1Brouha2003Merged$L1SeqFlank3p2x)
# 
# # Loop through flanking sequences above minimal length, locate them on the 
# # reference genome and get insert location descriptors
# for (i in which(FlankTotal > 100 & is.na(L1Brouha2003Merged$start_HG38))){
#   cat("Analyzing insertion", i, "of", nrow(L1Brouha2003Merged), "\n")
#   if (L1Brouha2003Merged$Strand[i] == "+") {
#     LPattern <- L1Brouha2003Merged$L1SeqFlank5p2x[i]
#     RPattern <- L1Brouha2003Merged$L1SeqFlank3p2x[i]
#     InsertionSite <- paste(L1Brouha2003Merged$L1SeqFlank5p[i],
#                            L1Brouha2003Merged$L1SeqFlank3p[i], sep = "")
#   } else {
#     DNASt5P   <- DNAString(L1Brouha2003Merged$L1SeqFlank5p[i])
#     DNASt3P   <- DNAString(L1Brouha2003Merged$L1SeqFlank3p[i])
#     DNASt5P2x <- DNAString(L1Brouha2003Merged$L1SeqFlank5p2x[i])
#     DNASt3P2x <- DNAString(L1Brouha2003Merged$L1SeqFlank3p2x[i])
#     LPattern  <- reverseComplement(DNASt3P2x) 
#     RPattern  <- reverseComplement(DNASt5P2x) 
#     InsertionSite <- paste(reverseComplement(DNASt3P),
#                            reverseComplement(DNASt5P), sep = "")
#   }
#   Chrom    <- paste(L1Brouha2003Merged$Chr[i])
#   ChromSeq <- BSgenome.Hsapiens.UCSC.hg38[[Chrom]]
#   String   <- matchLRPatterns(LPattern, RPattern, max.gaplength = 400, ChromSeq,
#                               max.Lmismatch = 10, max.Rmismatch = 10,
#                               with.Lindels = T, with.Rindels = T)
#   InsertStart <- NA
#   if (length(String) > 0){
#     InsertDF$FlankStart[i]        <- start(String) 
#     InsertDF$InsertionStartRel[i] <- 200 
#     InsertDF$InsertionStartAbs[i] <- start(String) + 200
#     InsertDF$NrNuc5p[i]           <- 0 
#     InsertDF$NrNuc3p[i]           <- 0 
#     pwA <- pairwiseAlignment(InsertionSite, String, type = "local")
#     Indel <- subsetByOverlaps(unlist(pwA@subject@indel), IRanges(start = 100, 
#                                                                  end = 100))
#     if (length(Indel) > 0){
#       InsertDF$InsertionStartRel[i] <- 100 + start(Indel) - 1 
#       InsertDF$InsertionStartAbs[i] <- start(String) + 100 + start(Indel) - 1
#       InsertDF$NrNuc5p[i]           <- 100 - start(Indel)
#       InsertDF$NrNuc3p[i]           <- end(Indel) - 100
#       
#     }
#   }
# }
# 
# # Re-run previous analysis for all insertions that could not be matched to 
# # reference but this time allowing for indels in flanking sequences
# for (i in which(FlankTotal > 100 & is.na(InsertDF$InsertionStartRel)) ){
#   cat("Analyzing insertion", i, "of", nrow(Beck2010TableWithL1), "\n")
#   if (Beck2010TableWithL1$Strand[i] == "+") {
#     LPattern <- Beck2010TableWithL1$L1SeqFlank5p2x[i]
#     RPattern <- Beck2010TableWithL1$L1SeqFlank3p2x[i]
#     InsertionSite <- paste(Beck2010TableWithL1$L1SeqFlank5p[i],
#                            Beck2010TableWithL1$L1SeqFlank3p[i], sep = "")
#   } else {
#     DNASt5P   <- DNAString(Beck2010TableWithL1$L1SeqFlank5p[i])
#     DNASt3P   <- DNAString(Beck2010TableWithL1$L1SeqFlank3p[i])
#     DNASt5P2x <- DNAString(Beck2010TableWithL1$L1SeqFlank5p2x[i])
#     DNASt3P2x <- DNAString(Beck2010TableWithL1$L1SeqFlank3p2x[i])
#     LPattern  <- reverseComplement(DNASt3P2x) 
#     RPattern  <- reverseComplement(DNASt5P2x) 
#     InsertionSite <- paste(reverseComplement(DNASt3P),
#                            reverseComplement(DNASt5P), sep = "")
#   }
#   Chrom    <- paste("chr", Beck2010TableWithL1$Chromosome[i], sep = "")
#   ChromSeq <- BSgenome.Hsapiens.UCSC.hg38[[Chrom]]
#   String   <- matchLRPatterns(LPattern, RPattern, max.gaplength = 300, ChromSeq,
#                               max.Lmismatch = 5, max.Rmismatch = 5, 
#                               with.Lindels = T, with.Rindels = T)
#   InsertStart <- NA
#   if (length(String) > 0){
#     InsertDF$FlankStart[i]        <- start(String) 
#     InsertDF$InsertionStartRel[i] <- 200 
#     InsertDF$InsertionStartAbs[i] <- start(String) + 200
#     InsertDF$NrNuc5p[i]           <- 0 
#     InsertDF$NrNuc3p[i]           <- 0 
#     pwA <- pairwiseAlignment(InsertionSite, String, type = "local")
#     Indel <- subsetByOverlaps(unlist(pwA@subject@indel), IRanges(start = 100, 
#                                                                  end = 100))
#     if (length(Indel) > 0){
#       InsertStart                   <- start(String) + 100 + start(Indel) - 1
#       InsertDF$InsertionStartRel[i] <- 100 + start(Indel) - 1 
#       InsertDF$InsertionStartAbs[i] <- start(String) + 200
#       InsertDF$NrNuc5p[i]           <- 100 - start(Indel)
#       InsertDF$NrNuc3p[i]           <- end(Indel) - 100
#       
#     }
#   }
# }

sum(is.na(InsertDF$InsertionStartAbs))
InsertDFBeck2010 <- InsertDF
Beck2010TableWithL1InsertLoc <- cbind(Beck2010TableWithL1, InsertDF)
write.csv(Beck2010TableWithL1InsertLoc, "D:/L1polymORF/Data/Beck2010_mergedTable_withL1InsertLoc.csv")
Beck2010TableWithL1InsertLoc <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withL1InsertLoc.csv")

#################################################
#                                               #
#   Assemble L1 table from Seleme et al. 2006   #
#                                               #
#################################################

# Read in tables for individual loci
Seleme2006L1A <- read.csv("D:/L1polymORF/Data/Seleme2006L1A.csv", as.is = T)
Seleme2006L1B <- read.csv("D:/L1polymORF/Data/Seleme2006L1B.csv", as.is = T)
Seleme2006L1C <- read.csv("D:/L1polymORF/Data/Seleme2006L1C.csv", as.is = T)

# Get indices for each of the accession numbers
idxA <- which(L1Brouha2003Merged$SeqSource == Seleme2006L1A$AccessionNr[1])
idxB <- which(L1Brouha2003Merged$SeqSource == Seleme2006L1B$AccessionNr[1])
idxC <- which(L1Brouha2003Merged$SeqSource == Seleme2006L1C$AccessionNr[1])

# Create genomic ranges for hot L1
GRanges_HotL1 <- GRanges(seqnames = L1Brouha2003Merged$Chr[1:3], 
                         ranges = IRanges(start = L1Brouha2003Merged$start_HG38[1:3],
                                          end = L1Brouha2003Merged$end_HG38[1:3]),
                         strand = L1Brouha2003Merged$strand[1:3])
Seqs_HotL1 <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges_HotL1)


# Determine sequence match
Table <- Seleme2006L1B
idx    <- which(L1Brouha2003Merged$SeqSource == Table$AccessionNr[1])
Seq    <- as.character(Seqs_HotL1[idx][[1]])

ColNameX   <- substr(colnames(Table), 2, nchar(colnames(Table)))
ColNameNum <- as.numeric(ColNameX)
idxPolyPos <- which(!is.na(ColNameNum) & Table[1,] %in% c("A", "G", "C", "T"))
PolyPos    <- ColNameNum[idxPolyPos]
names(PolyPos) <- idxPolyPos
Table[1,idxPolyPos]

# Test that Sequence of Brouha and Seleme match
OffSetVals <- -100:100
PMatch <- sapply(OffSetVals, function(OffSet){
  SeqMatch <- sapply(idxPolyPos, function(x) {
    Table[1,x] == substr(Seq, PolyPos[as.character(x)] + OffSet, PolyPos[as.character(x)] + OffSet)
  })
  sum(SeqMatch) /length(SeqMatch)
  
})
plot(PMatch)
Table[1,idxPolyPos]
BestOff <- OffSetVals[which.max(PMatch)]  
max(PMatch)
x <- 5
SMatch <- sapply(idxPolyPos, function(x) {
  Table[1,x] == substr(Seq, PolyPos[as.character(x)] + BestOff, PolyPos[as.character(x)] + BestOff)
})
PolyPos[!SMatch]

# Define function to exctract sequence from Seleme Table
Seq   <- L1Brouha2003Merged$Seq[L1Brouha2003Merged$SeqSource == Acc]

# Test that Sequence of Brouha and Seleme match
PMatch <- sapply(OffSetVals, function(OffSet){
  SeqMatch <- sapply(idxPolyPos, function(x) {
    Table[1,x] == substr(Seq, PolyPos[as.character(x)] + OffSet, PolyPos[as.character(x)] + OffSet)
  })
  sum(SeqMatch) /length(SeqMatch)
  
})
max(PMatch)


# Get the sequence from

#################################################
#                                               #
#        Put together data from                 #
#        Brouha, Beck and Seleme                #
#                                               #
#################################################

# Replace start and end for non-reference insertion in Beck2010 data
Beck2010TableWithL1InsertLoc$start_HG38 <-
  Beck2010TableWithL1InsertLoc$InsertionStartAbs
Beck2010TableWithL1InsertLoc$end_HG38 <-
  Beck2010TableWithL1InsertLoc$start_HG38 + 1

# Replace start and end for non-reference insertion in Brouha2003 data
blnBrouhaNonRef <- is.na(L1Brouha2003Merged$start_HG38)
L1Brouha2003Merged$start_HG38[blnBrouhaNonRef] <- 
  InsertDFBrouha2003$InsertionStartAbs[blnBrouhaNonRef]
L1Brouha2003Merged$Strand[is.na(L1Brouha2003Merged$Strand)] <- 
  L1Brouha2003Merged$strand[is.na(L1Brouha2003Merged$Strand)]

# Rename column so that they are consistent between datasets
colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Allele_frequency."] <- "Allele_frequency"
colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Chr"] <- "Chromosome"
colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Act_L1rp"] <- "Activity"

# Add additional columns
L1Brouha2003Merged$Coriell_ID <- NA
Beck2010TableWithL1InsertLoc$Allele_frequency <- NA
Beck2010TableWithL1InsertLoc$Chromosome <- paste("chr",
    Beck2010TableWithL1InsertLoc$Chromosome, sep = "")
Beck2010TableWithL1InsertLoc$Reference <- "Beck2010"

# Create tables for merging and merge them
CommonCols <- c("Accession", "Chromosome", "Activity", "Allele_frequency", 
                "Reference", "Coriell_ID", "Strand", "start_HG38", "end_HG38", "L1Seq",  
                "L1SeqFlank5p", "L1SeqFlank3p", "L1SeqFlank5p2x",
                "L1SeqFlank3p2x")
BeckForMeging   <- Beck2010TableWithL1InsertLoc[, CommonCols]
#BeckForMeging$Chromosome <- substr(BeckForMeging$Chromosome, 4, nchar(BeckForMeging$Chromosome))
BrouhaForMeging <- L1Brouha2003Merged[, CommonCols] 

# Merge both data sets
L1Catalogue <- rbind(BrouhaForMeging, BeckForMeging)

# Write cataloge
TStamp <- gsub(" ", "_", date())
TStamp <- gsub(":", "-", TStamp)
CataloguePath <- paste("D:/L1polymORF/Data/L1Catalogue_", TStamp, ".csv", sep = "")
write.csv(L1Catalogue, CataloguePath)

#################################################
#                                               #
#         Find variable sites                   #
#                                               #
#################################################

if (blnFindVarSites){
  
  # Loop through sequences and align them to consensus sequence
  AlignList <- sapply(L1Brouha2003Merged$Seq, function(x){
    if (is.na(x)){
      NULL
    } else {
      pairwiseAlignment(x, L1Consens)
    }
  })
  
  # Get all mismatch positions
  MisMatchList <- lapply(AlignList, function(x) {
    if (is.null(x)){
      NULL
    } else {
      x@subject@mismatch[[1]]
    }
  })
  AllMisMatches <- unique(unlist(MisMatchList))
  
  # Summarize mismatches
  blnNotNULLSeq <- sapply(MisMatchList, function(x) !is.null(x))
  NrSNPsPerSeq <- sapply(MisMatchList[blnNotNULLSeq], length)
  hist(NrSNPsPerSeq, breaks = seq(0,500, 10),
       xlab = "Number SNPs per sequence")
  
  # Get a matrix that counts the number of mismatching sequences per position
  SeqChar     <- as.character(L1Brouha2003Merged$Seq)
  MisMatchMat <- matrix(nrow = length(MisMatchList), ncol = max(nchar(SeqChar)))
  idxVect     <- 1:ncol(MisMatchMat)
  for (i in 1:nrow(MisMatchMat)){
    MisMatchMat[i,] <- idxVect %in% MisMatchList[[i]]
  }
  
  NrMM <- colSums(MisMatchMat)
  plot(which(NrMM > 0), NrMM[NrMM > 0], xlab = "Position within L1",
       ylab = "Number of sequences with mismatch")
  MMsm <- supsmu(1:length(NrMM), 1*(NrMM > 0), span = 0.1)
  lines(MMsm$x, 100* MMsm$y, col = "red")
  
}

#################################################
#                                               #
#         Find insertion sites                  #
#         of non-reference L1                   #
#                                               #
#################################################



###############################################
#                                             #
#     Write out sequences and align them      #
#                                             #
###############################################

if (blnAlignsequences) {
  # Read prototypic L1 sequences (from http://www.girinst.org/repbase/)  
  L1seqs <- read.dna("D:/L1polymORF/Data/Homo_sapiens_L1", format = "fasta",
                     as.character = T)
  
  # Remove \t from names of prototypic L1 sequences
  names(L1seqs) <- gsub("\t", "_", names(L1seqs))
  names(L1seqs) <- gsub(" ", "_", names(L1seqs))
  SeqLengths    <- sapply(L1seqs, length)
  
  # Find L1PA with maximum length
  idxL1PA     <- grep("L1PA", names(L1seqs))
  idxMaxL1PA  <- idxL1PA[which.max(SeqLengths[idxL1PA])]
  
  # Add L1PA as root (first sequence)
  idxWithSeq <- !is.na(L1Brouha2003Merged$Seq)
  L1HS_SeqList <- lapply(L1Brouha2003Merged$Seq[idxWithSeq],
                         function(x) tolower(s2c(x)))
  L1HS_withRoot <- c(L1seqs[idxMaxL1PA], L1HS_SeqList)
  names(L1HS_withRoot) <- c("L1PA", L1Brouha2003Merged$SeqSource[idxWithSeq])
  sapply(L1HS_withRoot, length)
  # Write sequences as fasta file and align
  write.fasta(L1HS_withRoot, names(L1HS_withRoot), 
              file.out = "D:/L1polymORF/Data/L1SequencesBrouha2003Unaligned.fas")
  run_MUSCLE(InputPath = "D:/L1polymORF/Data/L1SequencesBrouha2003Unaligned.fas", 
             OutputPath = "D:/L1polymORF/Data/L1SequencesBrouha2003Aligned.fas")
  
  # Read alignment and write it out as nexus file
  L1HSAligned <- read.dna("D:/L1polymORF/Data/L1SequencesBrouha2003Aligned.fas",
                          format = "fasta", as.character = T)
  write.nexus.data(L1HSAligned, "D:/L1polymORF/Data/L1SequencesBrouha2003Aligned.nex")
  
  # Add zero for the activitiy of L1PA
  Act_withRoot <- c(0, as.numeric(as.character(L1Brouha2003Merged$Act_L1rp[idxWithSeq])))
  names(Act_withRoot) <- names(L1HS_withRoot)
  
  # Code values in three classes
  Act_BayesTraits <- Act_withRoot
  Act_BayesTraits[Act_withRoot == 0]                   <- 0L
  Act_BayesTraits[Act_withRoot > 0 & Act_withRoot < 5] <- 1L
  Act_BayesTraits[Act_withRoot >= 5]                   <- 2L
  
  # Save Bayestraits file as text
  write.table(Act_BayesTraits, file = "D:/L1polymORF/Data/ActivityBT_Brouha2003.txt",
              col.names = F, quote = F)
  
}


