##############################################
#
# General description:
#
#   The following script reads in data sources on full-length L1 insertions and
#   locates the insertions in the reference genome

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
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Length of sequences flanking L1 insertions [bp] used to locate insertions
FlankSize <- 100

# Length of an L1 insertion with flanking sequences
LengthL1WithFlank <- 20000

# Set of columns in L1 catalog (column names are used throughout code so any 
# changes in column names have to be replaced throughout the code)
CommonCols <- c(
  "Accession",         # Accession number of clone sequence *
  "Allele",            # Allele ID (> 1 for multiple alleles of same insertion)
  "Chromosome",        # Chromosome of L1 insertion
  "Activity",          # Measured insertion activity (%)
  "Allele_frequency",  # Allele frequency of insertion
  "Reference",         # Paper publishing insertion
  "Coriell_ID", 
  "Strand",            # Standedness of L1 relative to clone sequence
  "start_HG38",        # start of L1 within hg38
  "end_HG38",          # end of L1 within hg38
  "L1SeqSourceType",   # Type of source for L1 sequence (BAC vs Fosmid)
  "L1Seq",             # L1 sequence
  "L1SeqFlank5p", 
  "L1SeqFlank3p", 
  "L1SeqFlank5p2x",
  "L1SeqFlank3p2x", 
  "start_Clone",       # start of L1 within clone sequence
  "end_Clone",         # end of L1 within clone sequence
  "strand_ClonetoRef", # Standedness of clone relative to reference sequence
  "strand_L1toRef"     # Standedness of L1 relative to reference sequence
)

# Boolean indicators for whether to perform particular processes
blnBuildBrouha2003 <- T
blnBuildBeck2010   <- T
blnBuildSeleme2006 <- T
blnMergeTables     <- T
blnFindVarSites    <- F
blnAlignsequences  <- F

#######################################
#                                     #
#   Define auxilliary functions       #
#                                     #
#######################################

# Function to create empty columns to append to an existing data table
CreateEmptyCols <- function(DataTable){
  ColsToAppend <- setdiff(CommonCols, colnames(DataTable))
  AppendData   <- matrix(nrow = nrow(DataTable), ncol = length(ColsToAppend))
  AppendData   <- as.data.frame(AppendData)
  colnames(AppendData) <- ColsToAppend
  AppendData
}

# Define function to switch strand
StrandSwitch <- function(Strand) switch(Strand, '+' = "-", '-' = "+")

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
GRanges_L1repMask_Hg38_Large <- GenomicRanges::resize(GRanges_L1repMask_Hg38, 
                                       width = LengthL1WithFlank, fix = "center")

# Get the 500 nuc upstream and downstream of the L1
FlankLeft <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                     names  = seqnames(GRanges_L1repMask_Hg38), 
                     start  = start(GRanges_L1repMask_Hg38)  - FlankSize, 
                     end    = start(GRanges_L1repMask_Hg38)) 
FlankRight <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                     names  = seqnames(GRanges_L1repMask_Hg38), 
                     start  = end(GRanges_L1repMask_Hg38), 
                     end    = end(GRanges_L1repMask_Hg38) + FlankSize) 
FlankLeft_RC  <- reverseComplement(FlankLeft)
FlankRight_RC <- reverseComplement(FlankRight)


L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges_L1repMask_Hg38)

#################################################
#                                               #
#   Process L1 table form Brouha et al. 2003    #
#                                               #
#################################################

if (blnBuildBrouha2003){
  
  # Read in table of hot L1 (obtained from Brouha et al.2003 PNAS)
  cat("******  Process L1 table by Brouha et al 2003  ******* \n\n")
  L1Brouha2003_raw <- read.delim('D:/L1polymORF/Data/L1HotTable_raw.txt', 
                                 sep = " ", as.is = T)

  # Find duplicated column names 
  CNames        <- colnames(L1Brouha2003_raw)
  ColNameSuffix <- substr(CNames, nchar(CNames) - 1, nchar(CNames))
  blnDuplCols   <- ColNameSuffix == ".1"
  
  # Separate table according to duplicated column names and append duplicated part
  L1Brouha2003Dupl <- L1Brouha2003_raw[,blnDuplCols]
  colnames(L1Brouha2003Dupl) <- CNames[!blnDuplCols]
  L1Brouha2003Table    <- rbind(L1Brouha2003_raw[,!blnDuplCols],L1Brouha2003Dupl)
  
  # Add chromosome column
  LChar          <- as.character(L1Brouha2003Table$Locus)
  LocusSpl1      <- sapply(LChar, function(x) strsplit(x, "[p-q]")[[1]][1])
  L1Brouha2003Table$Chr <- paste("chr", LocusSpl1, sep = "")
  
  # Rename columns
  colnames(L1Brouha2003Table)[colnames(L1Brouha2003Table) == "Accession_no."] <- "Accession"
  
  # Add missing empty columns
  ColsToAdd <- CreateEmptyCols(L1Brouha2003Table)
  L1Brouha2003Table <- cbind(L1Brouha2003Table, ColsToAdd)

  # Fill in data
  L1Brouha2003Table$Allele  <- 1
  L1Brouha2003Table$L1SeqSourceType  <- "BAC"
  L1Brouha2003Table$Reference        <- "Brouha2003"
  
  # Download sequences 
  cat("Get sequences for L1 loci by Brouha et al 2003\n")
  load("D:/L1polymORF/Data/Brouha2003Sequences.RData")

  # Re-write sequences in upper case and turn it into DNAStringset
  L1SeqHot <- toupper(L1SeqHot)
  L1SeqHotDNAStSet   <- DNAStringSet(L1SeqHot)
  L1SeqHotDNAStSetRV <- reverseComplement(L1SeqHotDNAStSet)
  
  # Get strand and chromosome from repeatMasker table
  Strand_L1repMask_Hg38 <- as.vector(strand(GRanges_L1repMask_Hg38))
  chr_L1repMask_Hg38    <- as.vector(seqnames(GRanges_L1repMask_Hg38))
  
  # Loop through flanking sequences and find parts of the BAC clones that match
  # flanking sequences
  for (i in 1:length(FlankLeft)){
    cat("Processing L1", i, "of", length(FlankLeft), "\n")
    Seq   <- L1HSSeq[i]

    # Get rows in Brouha table that match the current chromosome and subset the 
    # clone sequences
    idxChr   <- which(L1Brouha2003Table$Chr == chr_L1repMask_Hg38[i])
    DNAStSet <- L1SeqHotDNAStSet[idxChr]
    
    # Loop through all the clone sequences with matching chromosome and 
    # determine whether flanks of L1 in reference genome match clone sequence
    FlankMatches <- lapply(DNAStSet, function(x){
      matchLRPatterns(FlankLeft[[i]], FlankRight[[i]], max.gaplength = 6500, 
                      subject = x, max.Lmismatch = 5, max.Rmismatch = 5,
                      with.Lindels = T, with.Rindels = T)
    })
    idxMatch <- which(sapply(FlankMatches, length) > 0)
    CloneStrand <- "+"

    # Look for match of reverse complement if the original did not match
    if (length(idxMatch) == 0) {
      FlankMatches <- lapply(DNAStSet, function(x){
        matchLRPatterns(FlankRight_RC[[i]], FlankLeft_RC[[i]], max.gaplength = 6500,
                        subject = x, max.Lmismatch = 5, max.Rmismatch = 5,
                        with.Lindels = T, with.Rindels = T)
      })
      idxMatch <- which(sapply(FlankMatches, length) > 0)
      CloneStrand <- "-"
    }
    if (length(idxMatch) > 0) {
      mViews   <- FlankMatches[idxMatch][[1]]
      idxAbove6000 <- width(mViews) > 6000 & width(mViews) < 7000
      if (sum(idxAbove6000) == 1){
        S <- start(mViews[idxAbove6000])
        E <- end(mViews[idxAbove6000])
#        Seq   <- mViews[idxAbove6000]@subject[S:E] 
        Seq5P <- mViews[idxAbove6000]@subject[(S - FlankSize):(S - 1)] 
        Seq3P <- mViews[idxAbove6000]@subject[(E + 1):(E + FlankSize)] 
        L1Brouha2003Table$L1SeqFlank5p[idxChr[idxMatch]]  <- as.character(Seq5P)
        L1Brouha2003Table$L1SeqFlank3p[idxChr[idxMatch]]  <- as.character(Seq3P)
        L1Brouha2003Table$L1Seq[idxChr[idxMatch]]         <- as.character(Seq)
        L1Brouha2003Table$start_HG38[idxChr[idxMatch]]    <- start(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$end_HG38[idxChr[idxMatch]]      <- end(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$start_Clone[idxChr[idxMatch]]   <- S
        L1Brouha2003Table$end_Clone[idxChr[idxMatch]]     <- E
        L1Brouha2003Table$strand_ClonetoRef[idxChr[idxMatch]] <- CloneStrand
        L1Brouha2003Table$Strand[idxChr[idxMatch]]        <- Strand_L1repMask_Hg38[i]
      } else {
        cat("More than one set of flanking sequences enclose a stretch above 6000\n")
      }
      
    }
  }
  
  # Break L1Brouha2003Table into entries with and without L1 insertion on
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
                               MinMatchWidth = 5500, FlankSize = FlankSize,
                               Chromosomes = L1Brouha2003NonRef$Chr,
                               blnLocateL1inRef = T)

  # Replace each colum with the newly created dataframe
  L1DFBrouha2003$start_HG38 <- L1DFBrouha2003$start_Ref
  L1DFBrouha2003$end_HG38   <- L1DFBrouha2003$end_Ref
  for (Col in intersect(colnames(L1DFBrouha2003), colnames(L1Brouha2003NonRef)) ){
    L1Brouha2003NonRef[,Col] <- L1DFBrouha2003[,Col]
  }
  
  # Merge both tables
  L1Brouha2003Merged <- merge(L1Brouha2003Ref, L1Brouha2003NonRef, all = T)
  
  # Rename column so that they are consistent between datasets
  colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Allele_frequency."] <- "Allele_frequency"
  colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Chr"] <- "Chromosome"
  colnames(L1Brouha2003Merged)[colnames(L1Brouha2003Merged) == "Act_L1rp"] <- "Activity"
  
  # Write out table with sequence data
  cat("Saving table D:/L1polymORF/Data/L1Brouha2003.csv \n")
  write.csv(L1Brouha2003Merged, "D:/L1polymORF/Data/L1Brouha2003.csv",
            row.names = F)
  
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
  cat("******  Process L1 table by Beck et al 2010  ******* \n\n")
  Beck2010Table <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withEmptySite.csv",
                            as.is = T)
  Beck2010Table$Chromosome <- paste("chr", Beck2010Table$Chromosome, sep = "")
  Beck2010Table$Allele     <- 1
  L1DFBeck2010  <- getNonRefL1(L1Consens, AccNrs = Beck2010Table$Accession, 
                               MinMatchWidth = 5500, FlankSize = FlankSize,
                               Chromosomes = Beck2010Table$Chromosome,
                               blnLocateL1inRef = F)
  Beck2010TableWithL1 <- cbind(Beck2010Table, L1DFBeck2010)  
  Beck2010TableWithL1$start_HG38 <- Beck2010TableWithL1$start_Ref
  Beck2010TableWithL1$end_HG38   <- Beck2010TableWithL1$end_Ref
  cat("Saving table D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv \n")
  write.csv(Beck2010TableWithL1, "D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv",
            row.names = F)
} else {
  Beck2010TableWithL1 <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv",
                                  as.is = T)
}

#################################################
#                                               #
#   Assemble L1 table from Seleme et al. 2006   #
#                                               #
#################################################

if (blnBuildSeleme2006) {
  cat("******  Process L1 table by Seleme et al 2006  ******* \n\n")
  
  # Read in tables for individual loci
  Seleme2006L1A <- read.csv("D:/L1polymORF/Data/Seleme2006L1A.csv", as.is = T)
  Seleme2006L1B <- read.csv("D:/L1polymORF/Data/Seleme2006L1B.csv", as.is = T)
  Seleme2006L1C <- read.csv("D:/L1polymORF/Data/Seleme2006L1C.csv", as.is = T)
  Seleme2006L1B[1,which(Seleme2006L1B[1,] == "TRUE")] <- "T"
  
  # Function to get the best offset to match positions from Seleme 2006
  MatchSeqs <- function(SelemeTable, OffSetVals = -100:100) {
    
    # Find entry in Brouha table corresponding to Seleme table
    idx    <- which(L1Brouha2003Merged$Accession == SelemeTable$AccessionNr[1])
    Seq    <- s2c(L1Brouha2003Merged$L1Seq[idx])
    
    # Get all polymorphic positions
    ColNameX   <- substr(colnames(SelemeTable), 2, nchar(colnames(SelemeTable)))
    ColNameNum <- as.numeric(ColNameX)
    idxNotPoly <- is.na(ColNameNum)
    idxPolyPos <- which(!is.na(ColNameNum) & SelemeTable[1,] %in% c("A", "G", "C", "T"))
    PolyPos    <- ColNameNum[idxPolyPos]
    names(PolyPos) <- idxPolyPos
    
    # Find best offset value that maximizes match between Seleme and Brouha sequences
    PMatch <- sapply(OffSetVals, function(OffSet){
      SeqMatch <- SelemeTable[1,idxPolyPos] == Seq[PolyPos + OffSet]
      sum(SeqMatch) / length(SeqMatch)
    })
    plot(PMatch)
    BestOff <- OffSetVals[which.max(PMatch)]  
    
    # Determine which positions match 
    SMatch <- SelemeTable[1,idxPolyPos] == Seq[PolyPos + BestOff]
    
    # Create sequences for each allele
    Seqs <- sapply(1:nrow(SelemeTable), function(x){
      DiffPos <- which(SelemeTable[x,idxPolyPos] != "")
      Seq[PolyPos[DiffPos] + BestOff] <- SelemeTable[x,idxPolyPos[DiffPos]]
      paste(Seq, collapse = "")
    })

    # Return useful info
    return(list(BestOff = BestOff, NoMatchPos = PolyPos[!SMatch],
                NoMatchSeleme = SelemeTable[1,idxPolyPos[!SMatch]],
                NoMatchSeq = Seq[PolyPos[!SMatch] + BestOff],
                Seqs = Seqs, rowBrouha = idx, idxNotPoly = idxNotPoly,
                PropMatch = max(PMatch)))
  }
  
  # Get the best offset for each locus
  MatchListA <- MatchSeqs(Seleme2006L1A, OffSetVals = -100:100) 
  MatchListB <- MatchSeqs(Seleme2006L1B, OffSetVals = -100:100) 
  MatchListC <- MatchSeqs(Seleme2006L1C, OffSetVals = -100:100) 
  
  # Append column with L1 sequences and save resulting tables
  for (l in LETTERS[1:3]){
    SelTab <- eval(parse(text = paste("Seleme2006L1", l, sep = "")))
    MatchL <- eval(parse(text = paste("MatchList", l, sep = "")))
    SelTab$L1Seq <- MatchL$Seqs
    ColsToAppend <- setdiff(CommonCols, colnames(SelTab))
    AppendData   <- CreateEmptyCols(SelTab)
    for (i in 1:nrow(SelTab)){
      AppendData[i,] <- L1Brouha2003Merged[MatchL$rowBrouha, ColsToAppend]
    }
    SelTab <- cbind(SelTab, AppendData)
    assign(paste("Seleme2006L1extended", l, sep = ""), SelTab)
  }

  Seleme2006Combined <- rbind(
    Seleme2006L1extendedA[,-grep("X", colnames(Seleme2006L1extendedA))], 
    Seleme2006L1extendedB[,-grep("X", colnames(Seleme2006L1extendedB))], 
    Seleme2006L1extendedC[,-grep("X", colnames(Seleme2006L1extendedC))]) 
  Seleme2006Combined$Reference <- "Seleme2006"
  
  # Rename column so that they are consistent between datasets
  colnames(Seleme2006Combined)[colnames(Seleme2006Combined) == "AccessionNr"] <- "Accession"
  colnames(Seleme2006Combined)[colnames(Seleme2006Combined) == "Chr"]  <- "Chromosome"
  colnames(Seleme2006Combined)[colnames(Seleme2006Combined) == "Freq"] <- "Allele_frequency"
  
  cat("Saving table D:/L1polymORF/Data/Seleme2006Combined.csv \n")
  write.csv(Seleme2006Combined, "D:/L1polymORF/Data/Seleme2006Combined.csv",
            row.names = F)

} else {
  
  # Read in tables for individual loci
  Seleme2006Combined <- read.csv("D:/L1polymORF/Data/Seleme2006Combined.csv", 
                            as.is = T)
} 

#################################################
#                                               #
#        Put together data from                 #
#        Brouha, Beck and Seleme                #
#                                               #
#################################################

if(blnMergeTables){
  
  cat("******  Merge tables  ******* \n\n")

    # Add missing empty columns to Beck2010TableWithL1
  ColsToAdd <- CreateEmptyCols(Beck2010TableWithL1)
  Beck2010TableWithL1 <- cbind(Beck2010TableWithL1, ColsToAdd)
  
  # Fill in data
  Beck2010TableWithL1$L1SeqSourceType  <- "Fosmid"
  Beck2010TableWithL1$Reference        <- "Beck2010"
  

  # Create tables for merging and merge them
  BeckForMerging   <- Beck2010TableWithL1[, CommonCols]
  BrouhaForMerging <- L1Brouha2003Merged[, CommonCols] 
  SelemeForMerging <- Seleme2006Combined[, CommonCols] 
  SelemeForMerging <- SelemeForMerging[SelemeForMerging$Allele != 1, ]
  
  # Merge all data sets
  L1Catalogue <- rbind(BrouhaForMerging, BeckForMerging, SelemeForMerging)
  
  # Write cataloge
  TStamp <- gsub(" ", "_", date())
  TStamp <- gsub(":", "-", TStamp)
  CataloguePath <- paste("D:/L1polymORF/Data/L1Catalog_", TStamp, ".csv", sep = "")
  cat("Saving new catalog file", CataloguePath, "\n")
  write.csv(L1Catalogue, CataloguePath, row.names = F)
  
}

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
  SeqChar     <- as.character(L1Catalogue$Seq)
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
  idxWithSeq <- !is.na(L1Catalogue$L1Seq)
  L1HS_SeqList <- lapply(L1Catalogue$L1Seq[idxWithSeq],
                         function(x) tolower(s2c(x)))
  L1HS_withRoot <- c(L1seqs[idxMaxL1PA], L1HS_SeqList)
  names(L1HS_withRoot) <- c("L1PA", paste(L1Catalogue$Accession[idxWithSeq],
                                  L1Catalogue$Allele[idxWithSeq], sep = "_"))
  sapply(L1HS_withRoot, length)
  lapply(L1HS_withRoot, function(x) which(x == "I"))
  
  # Write sequences as fasta file and align
  write.fasta(L1HS_withRoot, names(L1HS_withRoot), 
              file.out = "D:/L1polymORF/Data/L1Catalogue_Unaligned_withRoot.fas")
  run_MUSCLE(InputPath = "D:/L1polymORF/Data/L1Catalogue_Unaligned_withRoot.fas", 
             OutputPath = "D:/L1polymORF/Data/L1Catalogue_Aligned_withRoot.fas")
  
  # Read alignment and write it out as nexus file
  L1HSAligned <- read.dna("D:/L1polymORF/Data/L1Catalogue_Aligned_withRoot.fas",
                          format = "fasta", as.character = T)
  write.nexus.data(L1HSAligned, "D:/L1polymORF/Data/L1Catalogue_Aligned_withRoot.nex")
  
  # Add zero for the activitiy of L1PA
  Act_withRoot <- c(0, as.numeric(as.character(L1Catalogue$Activity[idxWithSeq])))
  names(Act_withRoot) <- names(L1HS_withRoot)
  
  # Code values in three classes
  Act_BayesTraits <- Act_withRoot
  Act_BayesTraits[Act_withRoot == 0]                       <- 0L
  Act_BayesTraits[Act_withRoot > 0 & Act_withRoot < 50]    <- 1L
  Act_BayesTraits[Act_withRoot >= 50 & Act_withRoot < 100] <- 2L
  Act_BayesTraits[Act_withRoot >= 100]                     <- 3L
  
  # Save Bayestraits file as text
  write.table(Act_BayesTraits, file = "D:/L1polymORF/Data/ActivityBT_L1Catalogue.txt",
              col.names = F, quote = F)
  
}


