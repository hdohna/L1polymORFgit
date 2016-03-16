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

# Read euL1db tables
Fam_eul1db   <- read.delim("D:/HumanGenome/Data/Family_eul1db", skip = 5)
Ind_eul1db   <- read.delim("D:/HumanGenome/Data/Individuals_eul1db", skip = 5)
Met_eul1db   <- read.delim("D:/HumanGenome/Data/Methods_eul1db", skip = 5)
Stud_eul1db  <- read.delim("D:/HumanGenome/Data/Study_eul1db", skip = 5)
Sampl_eul1db <- read.delim("D:/HumanGenome/Data/Samples_eul1db", skip = 5)
MRIP_eul1db  <- read.delim("D:/HumanGenome/Data/MRIP_eul1db", skip = 5)
SRIP_eul1db  <- read.delim("D:/HumanGenome/Data/SRIP_eul1db", skip = 5)

# Numeric reference coordinates for SRIP
SRIP_eul1db$RefStart <- as.numeric(as.character(SRIP_eul1db$ref.start))
SRIP_eul1db$RefStop <- as.numeric(as.character(SRIP_eul1db$ref_stop))
pmin(SRIP_eul1db$RefStart, SRIP_eul1db$RefStop)

# Read table from the dbRIP database
dbRIP <- read.delim("D:/L1polymORF/Data/L1_hg19_v2h.txt", header = F)
dbRIP <- dbRIP[,-24]
colnames(dbRIP) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score",
                     "strand","originalId", "forwardPrimer", "reversePrimer", "polyClass", 
                     "polyFamily", "polySubfamily", "polySeq", "polySource", "reference", 
                     "ascertainingMethod", "remarks", "tm", "fillsize", "emptysize", "disease",
                     "genoRegion")

# Read in repeatMasker table 
L1repMask_Hg19   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv")
L1repMask_Hg38   <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table.csv")

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
GRangesFromTable <- function(DataTable, SeqNCol, StartCol, EndCol){
  NewStart <- pmin(DataTable[,StartCol], DataTable[,EndCol])
  NewEnd   <- pmax(DataTable[,StartCol], DataTable[,EndCol])
  idxNoNA  <- (!is.na(NewStart)) & (!is.na(NewEnd))
  IR <- IRanges(start = NewStart[idxNoNA], end = NewEnd[idxNoNA])
  GRanges(seqnames = DataTable[idxNoNA,SeqNCol], ranges = IR)
}

# Genomic ranges for each dataset
GRanges_SRIPg      <- GRangesFromTable(SRIP_eul1db, "chromosome", "g_start", "g_stop")
GRanges_SRIP       <- GRangesFromTable(SRIP_eul1db, "chromosome", "RefStart", "RefStop")
GRanges_MRIP      <- GRangesFromTable(MRIP_eul1db, "chromosome", "start", "stop")
GRanges_dbRIP     <- GRangesFromTable(dbRIP, "chrom", "chromStart", "chromEnd")
GRanges_L1repMask_Hg19 <- GRangesFromTable(L1repMask_Hg19, "genoName", "genoStart", 
                                           "genoEnd")
GRanges_L1repMask_Hg38 <- GRangesFromTable(L1repMask_Hg38, "genoName", "genoStart", 
                                      "genoEnd")
GRanges_IntactL1 <- GRangesFromTable(IntactL1, "Chr_hg18", "Start_hg18", "End_hg18")

# Get ranges in hg19
GRanges_IntactL1hg19 <- liftOver(GRanges_IntactL1, 
         chain = import.chain("D:/L1polymORF/Data/hg17ToHg19.over.chain"))

NrMapped <- sapply(GRanges_IntactL1hg19, length)
sum(NrMapped == 1)

# Retain only the coordinates that are uniquely mapped and dd coordinates to 
# table
blnMapped <- NrMapped == 1
GRanges_IntactL1hg19 <- unlist(GRanges_IntactL1hg19[NrMapped == 1])
IntactL1$Chr_hg19[blnMapped]   <- as.vector(seqnames(GRanges_IntactL1hg19))
IntactL1$Start_hg19[blnMapped] <- start(GRanges_IntactL1hg19)
IntactL1$End_hg19[blnMapped]   <- end(GRanges_IntactL1hg19)

# Check how many intactL1s are in each dataset
blnIntactL1_SRIPg  <- overlapsAny(GRanges_IntactL1hg19, GRanges_SRIPg)
sum(blnIntactL1_SRIPg)
blnIntactL1_SRIP  <- overlapsAny(GRanges_IntactL1hg19, GRanges_SRIP)
sum(blnIntactL1_SRIP)
blnIntactL1_MRIP  <- overlapsAny(GRanges_IntactL1hg19, GRanges_MRIP, minoverlap = 200)
sum(blnIntactL1_MRIP)
nearest(GRanges_IntactL1hg19, GRanges_MRIP)
GRanges_MRIP[nearest(GRanges_IntactL1hg19, GRanges_MRIP)]
blnIntactL1_dbRIP <- overlapsAny(GRanges_IntactL1hg19, GRanges_dbRIP)
sum(blnIntactL1_dbRIP)
sum(blnIntactL1_dbRIP & blnIntactL1_MRIP)
blnIntactL1_L1repMask_Hg19 <- overlapsAny(GRanges_IntactL1hg19, GRanges_L1repMask_Hg19, minoverlap = 6000)
sum(blnIntactL1_L1repMask_Hg19)
sum(blnIntactL1_L1repMask_Hg19 & blnIntactL1_MRIP)

# Check how many repMask_Hg19  are in each dataset
blnL1repMask_SRIPg  <- overlapsAny(GRanges_L1repMask_Hg19, GRanges_SRIPg)
sum(blnL1repMask_SRIPg)
blnL1repMask_SRIP  <- overlapsAny(GRanges_L1repMask_Hg19, GRanges_SRIP)
sum(blnL1repMask_SRIP)
blnL1repMask_MRIP  <- overlapsAny(GRanges_L1repMask_Hg19, GRanges_MRIP)
sum(blnL1repMask_MRIP)
blnL1repMask_dbRIP <- overlapsAny(GRanges_L1repMask_Hg19, GRanges_dbRIP)
sum(blnL1repMask_dbRIP)
blnL1repMask_L1repMask_Hg38 <- overlapsAny(GRanges_L1repMask_Hg19, GRanges_L1repMask_Hg38)
sum(blnL1repMask_L1repMask_Hg38)
