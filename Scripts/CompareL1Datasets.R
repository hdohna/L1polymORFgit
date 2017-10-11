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

###############################
#                             #
#      Read in data           #
#                             #
###############################

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
MRIP        <- read.delim("D:/L1polymORF/Data/eul1db_MRIP.txt", skip = 5)
SRIP        <- read.delim("D:/L1polymORF/Data/eul1db_SRIP.txt", skip = 5)

# Read in table by Kuhn et al.2014
L1Kuhn2014 <- read.csv("D:/L1polymORF/Data/Kuhn et al 2014 PNAS Table S2.csv")
L1Kuhn2014_GR <- makeGRangesFromDataFrame(L1Kuhn2014)

# Read in table by Ewing et al. 2011
L1Ewing2011    <- read.csv("D:/L1polymORF/Data/Ewing2011Table_S2.csv")
L1Ewing2011_S4 <- read.csv("D:/L1polymORF/Data/Ewing2011Table_S4.csv")
L1Ewing2011_S4$chromosome <- paste("chr", L1Ewing2011_S4$Chr, sep = "")
L1Ewing2011$chromosome <- paste("chr", L1Ewing2011$Chr, sep = "")
L1Ewing2011_S4_GR <- makeGRangesFromDataFrame(L1Ewing2011_S4, ignore.strand = T,
                                           seqnames.field = "chromosome",
                                           start.field = "Window.start",
                                           end.field = "Window.stop")
L1Ewing2011_GR <- makeGRangesFromDataFrame(L1Ewing2011, ignore.strand = T,
                                           seqnames.field = "chromosome",
                                           start.field = "Avg..Loc.",
                                           end.field = "Avg..Loc.")
L1Ewing2011_GR_hg19 <- UniqueLiftover(L1Ewing2011_GR,
   ChainFilePath = "D:/L1polymORF/Data/hg18ToHg19.over.chain")$LiftedRanges

grep("10", L1Ewing2011$Avg..Loc.)  
###############################
#                             #
#      Process data           #
#                             #
###############################

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

###############################
#                             #
#      Check overlap          #
#                             #
###############################

##############
# Check overlap with MRIP
##############

# Find overlap between L1 catalog and MRIP
blnL1CatInMRIP <- overlapsAny(L1CatalogGR_hg19, MRIP_GR)
cat("Proportion of catalog L1 in MRIP:", 
    sum(blnL1CatInMRIP) /length(blnL1CatInMRIP), "\n")
blnCatInRef <- L1CatalogL1Mapped$blnInRef
table(blnL1CatInMRIP, blnCatInRef)
chisq.test(blnL1CatInMRIP, blnCatInRef)

# Find overlap between 1000 genome L1  and MRIP
blnL1_100GInMRIP <- overlapsAny(L1_1000G_GR_hg19, MRIP_GR)
cat("Proportion of 1000 Genome L1 in MRIP:", 
    sum(blnL1_100GInMRIP) /length(blnL1_100GInMRIP), "\n")
cat("Proportion of MRIP L1 in 1000 Genomex :", 
    sum(blnL1_100GInMRIP) /length(MRIP_GR), "\n")
length(L1_1000G_GR_hg19)
L1_1000G_reduced$InsLength[blnL1_100GInMRIP]

##############
# Check overlap with SRIP
##############

# Find overlap between L1 catalog and SRIP
blnL1CatInSRIP <- overlapsAny(L1CatalogGR_hg19, SRIP_GR)
cat("Proportion of catalog L1 in SRIP:", 
    sum(blnL1CatInSRIP) /length(blnL1CatInSRIP), "\n")
blnCatInRef <- L1CatalogL1Mapped$blnInRef
table(blnL1CatInSRIP, blnCatInRef)
chisq.test(blnL1CatInSRIP, blnCatInRef)

# Find overlap between 1000 genome L1  and SRIP
blnL1_100GInSRIP <- overlapsAny(L1_1000G_GR_hg19, SRIP_GR)
cat("Proportion of 1000 Genome L1 in SRIP:", 
    sum(blnL1_100GInSRIP) /length(blnL1_100GInSRIP), "\n")
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

##############
# Check overlap of L1base with others
##############

blnL1CatInL1base <- overlapsAny(L1CatalogGR, L1base_GR)
length(L1base_GR)
length(L1CatalogGR)
L1base_GRList_hg19 <- liftOver(L1base_GR, 
                           import.chain("D:/L1polymORF/Data/hg18ToHg19.over.chain"))
NrMapped <- sapply(L1base_GRList_hg19, length)
idxUnique <- which(NrMapped == 1)
L1base_GR_hg19   <- unlist(L1base_GRList_hg19[idxUnique])
blnL1baseIn1000G <- overlapsAny(L1base_GR_hg19, L1_1000G_GR_hg19)

##############
# Check overlap of L1Kuhn2014_GR with 1000 Genomes and MRIP
##############

blnL1_KuhnIn1000G <- overlapsAny(L1Kuhn2014_GR, L1_1000G_GR_hg19)
sum(blnL1_KuhnIn1000G)
blnL1_KuhnInMRIP <- overlapsAny(L1Kuhn2014_GR, MRIP_GR)
sum(blnL1_KuhnInMRIP)
table(blnL1_KuhnIn1000G, blnL1_KuhnInMRIP)
sum(blnL1_KuhnInMRIP | blnL1_KuhnIn1000G)
table(L1Kuhn2014$KR_KNR)
hist(width(L1Kuhn2014_GR))

##############
# Check overlap of L1Ewing2011_GR with 1000 Genomes and MRIP
##############
L1Ewing2011_GR_Large <- resize(L1Ewing2011_GR_hg19, fix = "center",width = 500)
blnL1_EwingIn1000G <- overlapsAny(L1Ewing2011_GR_Large, L1_1000G_GR_hg19)
sum(blnL1_EwingIn1000G)
blnL1_EwingInMRIP <- overlapsAny(L1Ewing2011_GR_Large, MRIP_GR)
sum(blnL1_EwingInMRIP)
table(blnL1_EwingIn1000G, blnL1_EwingInMRIP)
sum(blnL1_EwingInMRIP | blnL1_EwingIn1000G)
table(L1Ewing2011$KR_KNR)
hist(width(L1Ewing2011_GR))

###############################
#                             #
#   Additional analysis       #
#                             #
###############################

# Add a column for the presence in 1000 genome data
L1Ewing2011$in1000G <- overlapsAny(L1Ewing2011_GR_Large, L1_1000G_GR_hg19)

# Test whether presence depends on insertion length
LogReg <- glm(in1000G ~ Length, data  = L1Ewing2011, family = binomial(link = "logit"))
summary(LogReg)

# Check proportion in 1000 genomes for different fragment length
blnShortFragm <- L1Ewing2011$Length <= 500
blnMedFragm  <- L1Ewing2011$Length >= 2000 & L1Ewing2011$Length <= 5000
blnLongFragm <- L1Ewing2011$Length >= 6000 
sum(blnShortFragm)
sum(blnMedFragm)
sum(blnLongFragm)

mean(L1Ewing2011$X1000.Genomes.Data[blnShortFragm] > 0)
mean(L1Ewing2011$X1000.Genomes.Data[blnMedFragm] > 0)
mean(L1Ewing2011$X1000.Genomes.Data[blnLongFragm] > 0)

# Group insertion length in units of 500 bp
InsLengthGroups <- seq(0, 6750, 750)
InsLengthMid  <- InsLengthGroups[-length(InsLengthGroups)] + 250
L1Ewing2011$InsLengthGroup <- cut(L1Ewing2011$Length,
                                       breaks = InsLengthGroups)
In1000GMean <- aggregate(in1000G ~ InsLengthGroup, data = L1Ewing2011, FUN = mean)
In1000GMin <- aggregate(in1000G ~ InsLengthGroup, data = L1Ewing2011, FUN = min)
In1000GMax <- aggregate(in1000G ~ InsLengthGroup, data = L1Ewing2011, FUN = max)
In1000GSum <- aggregate(in1000G ~ InsLengthGroup, data = L1Ewing2011, FUN = sum)
In1000GVar  <- aggregate(in1000G ~ InsLengthGroup, data = L1Ewing2011, FUN = var)
plot(InsLengthMid, In1000GMean$in1000G, 
     xlab = "L1 insertion length [bp]",
     ylab = "Proportion in 1000 genome")
AddErrorBars(MidX = InsLengthMid, MidY = In1000GMean$in1000G, 
             ErrorRange = sqrt(In1000GVar$in1000G / In1000GSum$in1000G * 
                                 In1000GMean$in1000G),
             TipWidth = 50)
L1_1000G_reduced$InsLengthGroup <- cut(L1_1000G_reduced$InsLength,
                                       breaks = InsLengthGroups)
# Get mean frequency
L1_1000G_reduced$InsLengthGroup <- cut(L1_1000G_reduced$InsLength,
                                       breaks = InsLengthGroups)
FreqMean <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = mean)
FreqMin <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = min)
FreqMax <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = max)
FreqSum <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = sum)
FreqVar  <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = var)
plot(InsLengthMid, FreqMean$Frequency, ylim = c(0, 0.05), 
     xlab = "L1 insertion length [bp]",
     ylab = "Frequency")
AddErrorBars(MidX = InsLengthMid, MidY = FreqMean$Frequency, 
             ErrorRange = sqrt(FreqVar$Frequency / FreqSum$Frequency * 
                                 FreqMean$Frequency),
             TipWidth = 50)
CreateDisplayPdf("D:/L1polymORF/Figures/L1InsertionLengthVsFrequency1000Genomes.pdf")

cor.test(FreqMean$Frequency, In1000GMean$in1000G)
