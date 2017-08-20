##############################################
#
# General description:
#
#   The following script reads in retortansposn data of the euL1db data base
#   see http://eul1db.unice.fr/

# Input:
#
#    Family_eul1db: data on families
#    Individuals_eul1db: data on individuals
#    Methods_eul1db: data on methods
#    Study_eul1db: data on studies
#    Samples_eul1db: data on samples
#    MRIP_eul1db: data on methods
#    SRIP_eul1db: data on methods

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

# Read in data
Families    <- read.delim("D:/L1polymORF/Data/eul1db_Family.txt", skip = 5)
Individuals <- read.delim("D:/L1polymORF/Data/eul1db_Individuals.txt", skip = 5)
Methods     <- read.delim("D:/L1polymORF/Data/eul1db_Methods.txt", skip = 5)
Studies     <- read.delim("D:/L1polymORF/Data/eul1db_Study.txt", skip = 5)
Samples     <- read.delim("D:/L1polymORF/Data/eul1db_Samples.txt", skip = 5)
MRIP        <- read.delim("D:/L1polymORF/Data/eul1db_MRIP.txt", skip = 5)
SRIP        <- read.delim("D:/L1polymORF/Data/eul1db_SRIP.txt", skip = 5)
str(SRIP)
table(SRIP$integrity)
hist(SRIP$allele_frequency)

# Add SRIP info to MRIP info
IDmatch <- match(MRIP$X.mrip_accession.no, SRIP$mrip)
MRIP$integrity <- SRIP$integrity[IDmatch]
MRIP$subgroup  <- SRIP$sub_group[IDmatch]
table(MRIP$integrity)

# Turn MRIP into genomic ranges
MRIP_GR <- makeGRangesFromDataFrame(MRIP)
SRIP_GR <- GRanges(seqnames = paste("chr", SRIP$chromosome, sep = ""), 
                   ranges = IRanges(start = pmin(SRIP$g_start, SRIP$g_stop),
                                    end = pmax(SRIP$g_start, SRIP$g_stop)))
hist(width(MRIP_GR), breaks = seq(0, 50000, 100), xlim = c(0, 7000))

##############
# Check overlap with MRIP
##############

# Find overlap between L1 catalog and MRIP
blnL1CatInMRIP <- overlapsAny(L1CatalogGR_hg19, MRIP_GR)
cat("Proportion of catalog L1 in MRIP:", sum(blnL1CatInMRIP) /length(blnL1CatInMRIP))
blnCatInRef <- L1CatalogL1Mapped$blnInRef
table(blnL1CatInMRIP, blnCatInRef)
chisq.test(blnL1CatInMRIP, blnCatInRef)

# Find overlap between 1000 genome L1  and MRIP
blnL1_100GInMRIP <- overlapsAny(L1_1000G_GR_hg19, MRIP_GR)
cat("Proportion of 1000 Genome L1 in MRIP:", sum(blnL1_100GInMRIP) /length(blnL1_100GInMRIP))
length(L1_1000G_GR_hg19)
L1_1000G_reduced$InsLength[blnL1_100GInMRIP]

##############
# Check overlap with SRIP
##############

# Find overlap between L1 catalog and SRIP
blnL1CatInSRIP <- overlapsAny(L1CatalogGR_hg19, SRIP_GR)
cat("Proportion of catalog L1 in SRIP:", sum(blnL1CatInSRIP) /length(blnL1CatInSRIP))
blnCatInRef <- L1CatalogL1Mapped$blnInRef
table(blnL1CatInSRIP, blnCatInRef)
chisq.test(blnL1CatInSRIP, blnCatInRef)

# Find overlap between 1000 genome L1  and SRIP
blnL1_100GInSRIP <- overlapsAny(L1_1000G_GR_hg19, SRIP_GR)
cat("Proportion of 1000 Genome L1 in SRIP:", sum(blnL1_100GInSRIP) /length(blnL1_100GInSRIP))
length(L1_1000G_GR_hg19)
L1_1000G_reduced$InsLength[blnL1_100GInSRIP]
table(blnL1_100GInSRIP, blnL1_100GInMRIP)

# Match MRIP data to 1000 Genome data
Overlap_MRIP_1000G <- findOverlaps(L1_1000G_GR_hg19, MRIP_GR)
plot(L1_1000G_reduced$Frequency[Overlap_MRIP_1000G@from],
     MRIP$pseudoallelefreq[Overlap_MRIP_1000G@to])
table(MRIP$integrity[Overlap_MRIP_1000G@to])

boxplot(pseudoallelefreq ~ integrity, data = MRIP, subset = subgroup == "L1-Ta")
aggregate(pseudoallelefreq ~ integrity, data = MRIP, subset = subgroup == "L1-Ta", FUN = mean)
table(MRIP$subgroup)

hist(MRIP$pseudoallelefreq[MRIP$integrity == "full-length"],
     breaks = seq(0, 1, 0.01))
hist(MRIP$pseudoallelefreq[MRIP$integrity == "5prime-truncated"],
     breaks = seq(0, 1, 0.01))
