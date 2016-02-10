##############################################
#
# General description:
#
#   The following script reads a repeat table downloaded from the genome
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

# Files and folders
RepeatFile          <- "D:/L1polymORF/Data/repeatsHg38"
TabOutfileName_L1HS <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"
SeqOutfileName_L1HS <- "D:/L1polymORF/Data/L1Sequences_reptab.fas"
TabOutfileName_L1   <- "D:/L1polymORF/Data/L1_repeat_table.csv"
SeqOutfileName_L1   <- "D:/L1polymORF/Data/L1AllSequences_reptab.fas"

#######################################
#                                     #
#     Read bed file and               #
#   save fasta file with sequences    #
#                                     #
#######################################

cat("Read and process repeatMasker table \n")

# Read repeat table
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]

######
# L1HS only
######

# Subset to get only L1HS rows and save
L1HSTable <- RepeatTable[RepeatTable$repName == "L1HS",]
idxFullLenghtL1HS <- which(abs(L1HSTable$genoEnd - L1HSTable$genoStart) >= 6000)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1HSIRanges <- IRanges(start = L1HSTable$genoStart,
                     end = L1HSTable$genoEnd)
L1HSGRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1HSIRanges,
                     strand = L1HSTable$strand)
L1HSGRanges <- L1HSGRanges[width(L1HSGRanges) >= 6000]

# Get all L1 sequences  
cat("Get sequences for L1 element \n")
L1HSSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1HSGRanges, as.character = T)
L1HSSeqNames <- paste(as.vector(seqnames(L1HSGRanges)), start(L1HSGRanges), end(L1HSGRanges),
      strand(L1HSGRanges), sep = "_")

# Write L1 sequences as fasta file
L1HSList <- lapply(L1HSSeq, s2c)
write.fasta(L1HSList, L1HSSeqNames, SeqOutfileName_L1HS)

######
# All L1s
######

# Subset to get only L1 rows and save
L1Table <- RepeatTable[RepeatTable$repFamily == "L1",]
idxL1HS_inFullTable <- which(abs(L1Table$genoEnd - L1Table$genoStart) >= 6000 &
                                L1Table$repName == "L1HS")

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                       end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                       strand = L1Table$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000]

# Get all L1 sequences  
L1Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GRanges, as.character = T)
L1SeqNames <- paste(as.vector(seqnames(L1GRanges)), start(L1GRanges), end(L1GRanges),
                      strand(L1GRanges), sep = "_")

# Write L1 sequences as fasta file
L1List <- lapply(L1Seq, s2c)
write.fasta(L1List, L1SeqNames, SeqOutfileName_L1)

#######################################
#                                     #
#     Process hot L1 table            #
#                                     #
#######################################

# Read in table of hot L1 (obtained from Brouha et al.2003 PNAS)
cat("Read table with hot L1s \n")
L1Hot_raw <- read.delim('D:/L1polymORF/Data/L1HotTable_raw.txt', sep = " ", as.is = T)

# Find duplicated column names 
CNames        <- colnames(L1Hot_raw)
ColNameSuffix <- substr(CNames, nchar(CNames) - 1, nchar(CNames))
blnDuplCols   <- ColNameSuffix == ".1"

# Separate table according to duplicated column names and append duplucated part
L1HotDupl <- L1Hot_raw[,blnDuplCols]
colnames(L1HotDupl) <- CNames[!blnDuplCols]
L1HotTable    <- rbind(L1Hot_raw[,!blnDuplCols],L1HotDupl)

# Download sequences 
cat("Get sequences for hot L1 loci \n")
choosebank("genbank")
L1SeqHot <- sapply(L1HotTable$Accession_no., function(AccNr){
  x <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
  Seq <- getSequence(x$req)
  paste(Seq[[1]], collapse = "")
})
closebank()

# Match each L1HS to a locus 
cat("Match hot L1 loci to L1s from repeatMasker\n")
L1SeqHot <- toupper(L1SeqHot)
GrepMatch <- sapply(L1HSSeq, function(x) {
  idxMatch <- grep(substr(x, 1, 500), L1SeqHot, ignore.case = T)
  if (length(idxMatch) == 0) NA else idxMatch
})

# Match info from L1SeqHot table to L1HSTable
L1HSTable[,c("Group", "Allele_frequency.", "Act_L1rp")] <- NA
L1HSTable[idxFullLenghtL1HS, c("Group", "Allele_frequency.", "Act_L1rp")] <- 
  L1HotTable[GrepMatch,c("Group", "Allele_frequency.", "Act_L1rp")]

# Match info from L1SeqHot table to L1HSTable
L1Table[,c("Group", "Allele_frequency.", "Act_L1rp")] <- NA
L1Table[idxL1HS_inFullTable, c("Group", "Allele_frequency.", "Act_L1rp")] <- 
  L1HotTable[GrepMatch,c("Group", "Allele_frequency.", "Act_L1rp")]

# Write out L1HSTable
write.csv(L1HSTable, TabOutfileName_L1HS)

# Write out L1Table
write.csv(L1Table, TabOutfileName_L1)
