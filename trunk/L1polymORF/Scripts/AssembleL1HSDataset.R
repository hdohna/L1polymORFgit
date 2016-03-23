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

# Boolean indicators for whether to perform particular processes
blnBuildBrouha2003 <- F
blnBuildBeck2010   <- F
blnFindVarSites    <- F
blnAlignsequences  <- F

# Files and folders
RepeatFile          <- "D:/L1polymORF/Data/repeatsHg38"
TabOutfileName_L1HS <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"
SeqOutfileName_L1HS <- "D:/L1polymORF/Data/L1Sequences_reptab.fas"
SeqOutfileName_L1HS_length100 <- "D:/L1polymORF/Data/L1HSSequences_reptab_L100.fas"
TabOutfileName_L1   <- "D:/L1polymORF/Data/L1_repeat_table.csv"
SeqOutfileName_L1   <- "D:/L1polymORF/Data/L1AllSequences_reptab.fas"

# Files and folders
dbRIP_Path <- "D:/L1polymORF/Data/L1_hg19_v2h.txt"
L1Repeat_path   <- "D:/L1polymORF/Data/L1_repeat_table.csv"

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Load consensus sequence
L1Consens <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
L1Consens <- paste(L1Consens[[1]], collapse = "")
L1Consens <- DNAString(L1Consens)

# Read euL1db tables
Fam_eul1db   <- read.delim("D:/HumanGenome/Data/Family_eul1db", skip = 5)
Ind_eul1db   <- read.delim("D:/HumanGenome/Data/Individuals_eul1db", skip = 5)
Met_eul1db   <- read.delim("D:/HumanGenome/Data/Methods_eul1db", skip = 5)
Stud_eul1db  <- read.delim("D:/HumanGenome/Data/Study_eul1db", skip = 5)
Sampl_eul1db <- read.delim("D:/HumanGenome/Data/Samples_eul1db", skip = 5)
MRIP_eul1db  <- read.delim("D:/HumanGenome/Data/MRIP_eul1db", skip = 5)
#SRIP_eul1db  <- read.delim("D:/HumanGenome/Data/SRIP_eul1db", skip = 5)

# Numeric reference coordinates for SRIP
# SRIP_eul1db$RefStart <- as.numeric(as.character(SRIP_eul1db$ref.start))
# SRIP_eul1db$RefStop <- as.numeric(as.character(SRIP_eul1db$ref_stop))
# pmin(SRIP_eul1db$RefStart, SRIP_eul1db$RefStop)

# Read table from the dbRIP database
dbRIP <- read.delim("D:/L1polymORF/Data/L1_hg19_v2h.txt", header = F)
dbRIP <- dbRIP[,-24]
colnames(dbRIP) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score",
                     "strand","originalId", "forwardPrimer", "reversePrimer", 
                     "polyClass", "polyFamily", "polySubfamily", "polySeq", 
                     "polySource", "reference", "ascertainingMethod", "remarks", 
                     "tm", "fillsize", "emptysize", "disease",
                     "genoRegion")

# Remove non-L1HS
dbRIP <- dbRIP[!dbRIP$polySubfamily %in% c("L1PA2", "L1PA3", "L1M2_orf2", 
                                           "L1PA5", "L1PA4", "L1M2", "L1PB4", 
                                           "L1MA8"),]

# Read in repeatMasker table 
# L1repMask_Hg18 <- read.delim("D:/L1polymORF/Data/repeatsHg18")
# L1repMask_Hg18 <- L1repMask_Hg18[nchar(as.character(L1repMask_Hg18$genoName)) <= 5, ]

L1repMask_Hg19   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv")
L1repMask_Hg38   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table.csv")
table(L1repMask_Hg19$repName)

# Read in table of intact L1 
IntactL1       <- read.csv("D:/L1polymORF/Data/IntactL1.csv", as.is = T)
IntactL1$Chr <- paste("chr", IntactL1$Chr, sep = "")

# Rename start, end, and chromosme columns to indicate they refer to hg18
blnCoordCols <- colnames(IntactL1) %in% c("Chr", "Start", "End")
colnames(IntactL1)[blnCoordCols] <- paste(colnames(IntactL1)[blnCoordCols], 
                                          "hg18", sep = "_")

# Add coordinate columns for hg19
IntactL1$Chr_hg19   <- NA
IntactL1$Start_hg19 <- NA
IntactL1$End_hg19   <- NA

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
# GRanges_SRIPg      <- GRangesFromTable(SRIP_eul1db, "chromosome", "g_start", "g_stop")
# GRanges_SRIP       <- GRangesFromTable(SRIP_eul1db, "chromosome", "RefStart", "RefStop")
GRanges_MRIP       <- GRangesFromTable(MRIP_eul1db, "chromosome", "start", "stop")
GRanges_dbRIP      <- GRangesFromTable(dbRIP, "chrom", "chromStart", "chromEnd")
GRanges_L1repMask_Hg19 <- GRangesFromTable(L1repMask_Hg19, "genoName", "genoStart", 
                                           "genoEnd", "strand")
GRanges_L1repMask_Hg38 <- GRangesFromTable(L1repMask_Hg38, "genoName", "genoStart", 
                                      "genoEnd", "strand")
GRanges_IntactL1 <- GRangesFromTable(IntactL1, "Chr_hg18", "Start_hg18", "End_hg18")

# Get ranges in hg19
GRanges_IntactL1hg19 <- liftOver(GRanges_IntactL1, 
         chain = import.chain("D:/L1polymORF/Data/hg17ToHg19.over.chain"))

NrMapped <- sapply(GRanges_IntactL1hg19, length)

# Retain only the coordinates that are uniquely mapped and dd coordinates to 
# table
blnMapped <- NrMapped == 1
cat(sum(blnMapped), "of", length(GRanges_IntactL1), 
    "intact L1 could be mapped to hg19\n")  
GRanges_IntactL1hg19 <- unlist(GRanges_IntactL1hg19[NrMapped == 1])
IntactL1$Chr_hg19[blnMapped]   <- as.vector(seqnames(GRanges_IntactL1hg19))
IntactL1$Start_hg19[blnMapped] <- start(GRanges_IntactL1hg19)
IntactL1$End_hg19[blnMapped]   <- end(GRanges_IntactL1hg19)

#######################################
#                                     #
#     Compare datasets                #
#                                     #
#######################################

# Get names of datasets
GRangesNames <- c("GRanges_IntactL1hg19", "GRanges_MRIP",
                  "GRanges_dbRIP", "GRanges_L1repMask_Hg19")
GRangesNamesNames <- c(GRanges_IntactL1hg19 = "L1Base", GRanges_MRIP = "euL1db",
                  GRanges_dbRIP = "dbRIP", GRanges_L1repMask_Hg19 = "repeatMask")

# Remove strange genomic ranges
RemoveRanges <- function(Ranges, LowerW = 0, UpperW = Inf){
  Ranges[width(Ranges) >= LowerW & width(Ranges) <= UpperW]
}

# # Plot histogram of widths
# par(mfrow = c(2, 2))
# GeneralBreaks <- seq(0, 65000, 100)
# for (GRN in GRangesNames){
#   R <- eval(parse(text = GRN))
#   hist(width(R), xlab = "Genomic range width", main = GRangesNamesNames[GRN],
#        breaks = GeneralBreaks, xlim = c(min(width(R)), max(width(R))))
# }
# #CreateDisplayPdf("D:/L1polymORF/Figures/L1WidthPerDataSet_beforeAdj.pdf")
# 
# GRanges_MRIP  <- RemoveRanges(GRanges_MRIP, UpperW = 12000)
# GRanges_dbRIP <- RemoveRanges(GRanges_dbRIP, LowerW = 6000, UpperW = 12000)
# GRanges_L1repMask_Hg19 <- RemoveRanges(GRanges_L1repMask_Hg19, LowerW = 6000, 
#                                        UpperW = 12000)
# # Plot histogram of widths
# par(mfrow = c(2, 2))
# GeneralBreaks <- seq(0, 65000, 100)
# for (GRN in GRangesNames){
#   R <- eval(parse(text = GRN))
#   hist(width(R), xlab = "Genomic range width", main = GRangesNamesNames[GRN],
#        breaks = GeneralBreaks, xlim = c(min(width(R)), max(width(R))))
# }
#CreateDisplayPdf("D:/L1polymORF/Figures/L1WidthPerDataSet_afterAdj.pdf")

# Create matrix of the proportion of genomic ranges overlapping by at least one
# bp between different datasets
# PropAnyOverlap <- sapply(GRangesNames, function(x){
#   sapply(GRangesNames, function(y){
#     R1 <- eval(parse(text = x))
#     R2 <- eval(parse(text = y))
#     sum(overlapsAny(R1, R2)) / length(R1)
#   })
# })
# rownames(PropAnyOverlap) <- GRangesNamesNames[rownames(PropAnyOverlap)]
# colnames(PropAnyOverlap) <- GRangesNamesNames[colnames(PropAnyOverlap)]
# PropAnyOverlap <- round(PropAnyOverlap, 3)
# write.csv(PropAnyOverlap, file = "D:/L1polymORF/Data/DataSetComparePropAny.csv")
# 
# # Create matrix of the average number of nucleotides shared between overlapping
# # genomic regins between different datasets
# MeanPropBpOverlap <- sapply(GRangesNames, function(x){
#   sapply(GRangesNames, function(y){
#     R1 <- eval(parse(text = x))
#     R2 <- eval(parse(text = y))
#     Os <- findOverlaps(R1, R2)
#     qHits <- Os@queryHits
#     sHits <- Os@subjectHits
#     idxDupl   <- duplicated(qHits)
#     DuplQuery <- unique(qHits[idxDupl])
#     qHits <- qHits[!idxDupl]
#     sHits <- sHits[!idxDupl]
#     Props <- width(pintersect(R1[qHits], R2[sHits])) / width(R1[qHits])
#     for (i in DuplQuery){
#       idx <- i %in% Os@queryHits
#       NewProps <- width(pintersect(R1[Os@queryHits[idx]], R2[Os@subjectHits[idx]])) / 
#         width(R1[Os@queryHits[idx]])
#       Props[qHits == i] <- max(NewProps)
#     }
#     mean(Props)
#   })
# })
# rownames(MeanPropBpOverlap) <- GRangesNamesNames[rownames(MeanPropBpOverlap)]
# colnames(MeanPropBpOverlap) <- GRangesNamesNames[colnames(MeanPropBpOverlap)]
# MeanPropBpOverlap <- round(MeanPropBpOverlap, 3)
# write.csv(MeanPropBpOverlap, file = "D:/L1polymORF/Data/DataSetCompareMeanProp.csv")

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
FlankSize <- 100
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
  L1Brouha2003_raw <- read.delim('D:/L1polymORF/Data/L1Brouha2003Table_raw.txt', sep = " ", as.is = T)
  
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
  colnames(L1Brouha2003Table)[colnames(L1Brouha2003Table) == "Accession_no."] <- "SeqSource"
  
  # Add more columns
  L1Brouha2003Table$Seq         <- NA
  L1Brouha2003Table$SeqFlank5p  <- NA
  L1Brouha2003Table$SeqFlank3p  <- NA
  L1Brouha2003Table$SeqSourceType  <- "BAC"
  L1Brouha2003Table$start_HG38  <- NA
  L1Brouha2003Table$end_HG38    <- NA
  L1Brouha2003Table$strand      <- NA
  L1Brouha2003Table$Reference  <- "Brouha2003"
  
  # Download sequences 
  cat("Get sequences for hot L1 loci \n")
  load("D:/L1polymORF/Data/Brouha2003Sequences.RData")
  # choosebank("genbank")
  # L1SeqHot <- sapply(L1Brouha2003Table$SeqSource, function(AccNr){
  #   x <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
  #   Seq <- getSequence(x$req)
  #   paste(Seq[[1]], collapse = "")
  # })
  # closebank()
  # names(L1SeqHot)
  # save(list = "L1SeqHot", file = "D:/L1polymORF/Data/Brouha2003Sequences.RData")
  
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
      matchLRPatterns(LeftP, RightP, max.gaplength = 6100, 
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
        L1Brouha2003Table$SeqFlank5p[idxChr[idxMatch]]  <- as.character(Seq5P)
        L1Brouha2003Table$SeqFlank3p[idxChr[idxMatch]]  <- as.character(Seq3P)
        L1Brouha2003Table$Seq[idxChr[idxMatch]] <- as.character(Seq)
        L1Brouha2003Table$start_HG38[idxChr[idxMatch]] <- start(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$end_HG38[idxChr[idxMatch]]   <- end(GRanges_L1repMask_Hg38)[i]
        L1Brouha2003Table$strand[idxChr[idxMatch]] <- Strand_L1repMask_Hg38[i]
      } else {
        cat("More than one set of flanking sequences enclose a stretch above 6000\n")
      }
      
    }
  }
  
  
  # Loop through sequences that have no coordinates in the reference genome
  # and get the sequences through pairwise alignment
  for (i in c(4, which(is.na(L1Brouha2003Table$Seq)))){
    print(i)
    SubjectSeq <- L1SeqHotDNAStSet[i]
    lpA <- pairwiseAlignment(L1Consens, SubjectSeq, type = "local")
    if (width(lpA@subject) > 5500){
      Strand <- "+"
    } else {
      SubjectSeq <- L1SeqHotDNAStSetRV[i]
      lpA        <- pairwiseAlignment(L1Consens, SubjectSeq, type = "local")
      Strand     <- "-"
    }
    if (width(lpA@subject) > 5500){
      L1Brouha2003Table$Seq[i]    <- as.character(lpA@subject)
      L1Brouha2003Table$strand[i] <- Strand
      S <- start(lpA@subject@range)
      E <- end(lpA@subject@range)
      Seq5P <- SubjectSeq[[1]][max(1, (S - FlankSize)):(S - 1)] 
      Seq3P <- SubjectSeq[[1]][(E + 1):min(length(SubjectSeq[[1]]), (E + FlankSize))] 
      L1Brouha2003Table$SeqFlank5p[idxChr[idxMatch]]  <- as.character(Seq5P)
      L1Brouha2003Table$SeqFlank3p[idxChr[idxMatch]]  <- as.character(Seq3P)
    } else {
      cat("L1 consensus could not be matched to BAC clone", 
          L1Brouha2003Table$SeqSource[i], "\n")
    }
  }
  
  # Write out table with sequence data
  write.csv(L1Brouha2003Table, "D:/L1polymORF/Data/DataBrouha2003.csv")
  
} else {
  L1Brouha2003Table <- read.csv("D:/L1polymORF/Data/DataBrouha2003.csv", as.is = T)
}

#################################################
#                                               #
#   Process L1 table from Beck et al. 2003      #
#                                               #
#################################################

if (blnBuildBeck2010) {
  Beck2010Table <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withEmptySite.csv")
  AccChar       <- as.character(Beck2010Table$AccessionNumber)
  AccNrCleaned  <- sapply(AccChar,function(x) strsplit(x, "\\.")[[1]][1])
  L1DFBeck2010  <- getNonRefL1(L1Consens, AccNrs = AccNrCleaned, 
                               MinMatchWidth = 5500)
  Beck2010TableWithL1 <- cbind(Beck2010Table, L1DFBeck2010)  
  write.csv(Beck2010TableWithL1, "D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv")
  
} else {
  Beck2010TableWithL1 <- read.csv("D:/L1polymORF/Data/Beck2010_mergedTable_withL1.csv")
  
}

#################################################
#                                               #
#   Assemble L1 table from Seleme et al. 2003   #
#                                               #
#################################################



# Read in tables for individual loci
Seleme2006L1A <- read.csv("D:/L1polymORF/Data/Seleme2006L1A.csv", as.is = T)
Seleme2006L1B <- read.csv("D:/L1polymORF/Data/Seleme2006L1B.csv", as.is = T)
Seleme2006L1C <- read.csv("D:/L1polymORF/Data/Seleme2006L1C.csv", as.is = T)

# Define function to exctract sequence from Seleme Table
Table <- Seleme2006L1C
Acc   <- Table$AccessionNr[1]
Seq   <- L1Brouha2003Table$Seq[L1Brouha2003Table$SeqSource == Acc]

ColNameX   <- substr(colnames(Table), 2, nchar(colnames(Table)))
ColNameNum <- as.numeric(ColNameX)
idxPolyPos <- which(!is.na(ColNameNum) & Table[1,] != "")
PolyPos    <- ColNameNum[idxPolyPos]
names(PolyPos) <- idxPolyPos

# Test that Sequence of Brouha and Seleme match
OffSetVals <- -300:300
OffSet <- 99
x <- idxPolyPos[1]
PMatch <- sapply(OffSetVals, function(OffSet){
  SeqMatch <- sapply(idxPolyPos, function(x) {
    Table[1,x] == substr(Seq, PolyPos[as.character(x)] + OffSet, PolyPos[as.character(x)] + OffSet)
  })
  sum(SeqMatch) /length(SeqMatch)
  
})
plot(PMatch)

OffSetVals[which.max(PMatch)]  
max(PMatch)


Table <- Seleme2006L1C
Acc   <- Table$AccessionNr[1]
choosebank("genbank")
x <- query(listname = "L1", paste("AC=", Acc, sep = ""))
SeqBAC <- getSequence(x$req)
SeqBAC <- paste(SeqBAC[[1]], collapse = "")
SeqBAC <- toupper(SeqBAC)
closebank()
SeqBACDNASt <- DNAString(SeqBAC)
SeqBACDNAStRV <- reverseComplement(SeqBACDNASt)
SeqBAC_L1 <- SeqBACDNAStRV[113051:11906]
SeqBAC_L1RV <- reverseComplement(SeqBAC_L1)
PMatch <- sapply(OffSetVals, function(OffSet){
  SeqMatch <- sapply(idxPolyPos, function(x) {
    Pos <- PolyPos[as.character(x)] + OffSet
    Nuc <- SeqBAC_L1[Pos:Pos]
    Table[1,x] == as.character(Nuc)
  })
  sum(SeqMatch) /length(SeqMatch)
})
max(PMatch)
#################################################
#                                               #
#         Find variable sites                   #
#                                               #
#################################################

if (blnFindVarSites){
  # Loop through sequences and align them to consensus sequence
  AlignList <- sapply(L1Brouha2003Table$Seq, function(x){
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
  SeqChar     <- as.character(L1Brouha2003Table$Seq)
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
  idxWithSeq <- !is.na(L1Brouha2003Table$Seq)
  L1HS_SeqList <- lapply(L1Brouha2003Table$Seq[idxWithSeq],
                         function(x) tolower(s2c(x)))
  L1HS_withRoot <- c(L1seqs[idxMaxL1PA], L1HS_SeqList)
  names(L1HS_withRoot) <- c("L1PA", L1Brouha2003Table$SeqSource[idxWithSeq])
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
  Act_withRoot <- c(0, as.numeric(as.character(L1Brouha2003Table$Act_L1rp[idxWithSeq])))
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


