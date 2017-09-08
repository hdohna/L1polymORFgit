##############################################
#
# General description:
#
#   The following script takes data collected by Kuhn et al. (2014, PNAS)
#   and adds infromation from other databases such as eul1db 
#   (http://eul1db.unice.fr/) and 1000 genomes 

# Input:
#
#    Family_eul1db: data on families
#    Individuals_eul1db: data on individuals
#    Methods_eul1db: data on methods
#    Study_eul1db: data on studies
#    Samples_eul1db: data on samples
#    MRIP_eul1db: data on methods
#    SRIP_eul1db: data on methods
#    Kuhn2014PNAS_TableS2.csv: data collected by Kuhn et al. 2014

# Output:
#   
#    : ...

##############################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Load L1 catalog GenomicRanges
load("D:/L1polymORF/Data/L1CatalogGRanges.RData")
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

# Read in Charite L1 data
L1base    <- read.csv("D:/L1polymORF/Data/L1baseCharite.csv")
L1base_GR <- GRanges(seqnames = paste("chr", L1base$Chr, sep = ""),
    ranges = IRanges(start = pmin(L1base$Start, L1base$End),
                     end   = pmax(L1base$Start, L1base$End)))

# Read in eul1db data
Families    <- read.delim("D:/L1polymORF/Data/eul1db_Family.txt", skip = 5)
Individuals <- read.delim("D:/L1polymORF/Data/eul1db_Individuals.txt", skip = 5)
Methods     <- read.delim("D:/L1polymORF/Data/eul1db_Methods.txt", skip = 5)
Studies     <- read.delim("D:/L1polymORF/Data/eul1db_Study.txt", skip = 5)
MRIP        <- read.delim("D:/L1polymORF/Data/eul1db_MRIP.txt", skip = 5)
SRIP        <- read.delim("D:/L1polymORF/Data/eul1db_SRIP.txt", skip = 5)

# Add SRIP info to MRIP info
IDmatch <- match(MRIP$X.mrip_accession.no, SRIP$mrip)
MRIP$integrity <- SRIP$integrity[IDmatch]
MRIP$subgroup  <- SRIP$sub_group[IDmatch]

# Turn MRIP into genomic ranges
MRIP_GR <- makeGRangesFromDataFrame(MRIP)

# Read in table by Kuhn et al. 
L1Kuhn2014    <- read.csv("D:/L1polymORF/Data/Kuhn2014PNAS_TableS2.csv")
L1Kuhn2014_GR <- makeGRangesFromDataFrame(L1Kuhn2014)

# Read in table for reference genome
L1Ref    <- read.csv("D:/L1polymORF/Data/repeatsHg19_L1HS.csv")
L1Ref_GR <- makeGRangesFromDataFrame(L1Ref, 
                                     seqnames.field = "genoName",
                                     start.field="genoStart",
                                     end.field="genoEnd")

##############
# Add infor from 1000 Genomes and MRIP
##############

# Add a column for full length
L1Kuhn2014$blnFull <- NA

# Get overlap with MRIP and add info 
OverlapKuhn_MRIP <- findOverlaps(L1Kuhn2014_GR, MRIP_GR)
blnFullMRIP      <- MRIP$integrity == "full-length"
blnFullMRIP[MRIP$integrity == "unknown"] <- NA
L1Kuhn2014$blnFull[OverlapKuhn_MRIP@from] <- blnFullMRIP[OverlapKuhn_MRIP@to]
cat("L1 status added for", sum(!is.na(L1Kuhn2014$blnFull)), "entries\n")

# Get overlap with 1000 genome and add info
OverlapKuhn_1000G <- findOverlaps(L1Kuhn2014_GR, L1_1000G_GR_hg19)
blnFull1000G      <- L1_1000G_reduced$InsLength > 6000
L1Kuhn2014$blnFull[OverlapKuhn_1000G@from] <- blnFull1000G[OverlapKuhn_1000G@to]
cat("L1 status added for", sum(!is.na(L1Kuhn2014$blnFull)), "entries\n")
length(OverlapKuhn_1000G@to)

# Get overlap with reference genome and add info
OverlapKuhn_Ref <- findOverlaps(L1Kuhn2014_GR, L1Ref_GR)
blnFullRef      <- width(L1Ref_GR) > 6000
L1Kuhn2014$blnFull[OverlapKuhn_Ref@from] <- blnFull1000G[OverlapKuhn_Ref@to]
cat("L1 status added for", sum(!is.na(L1Kuhn2014$blnFull)), "entries\n")

# Add column with insertion width
L1Kuhn2014$InsWidth <- width(L1Kuhn2014_GR)
aggregate(InsWidth ~ KR_KNR, data = L1Kuhn2014, FUN = max)
aggregate(frequency ~ KR_KNR, data = L1Kuhn2014, FUN = mean)/20

hist(L1Kuhn2014$frequency[L1Kuhn2014$KR_KNR == "KR"])
