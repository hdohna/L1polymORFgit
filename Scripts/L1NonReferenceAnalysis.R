##############################################
#
# General description:
#
#   The following script analyzes non-reference L1 detected by our analysis 
#   and compares them from data form the euL1db data base
#   (see http://eul1db.unice.fr/)

# Input:
#
#    Family_eul1db: data on families
#    Individuals_eul1db: data on individuals
#    Methods_eul1db: data on methods
#    Study_eul1db: data on studies
#    Samples_eul1db: data on samples
#    MRIP_eul1db: data on methods
#    SRIP_eul1db: data on methods
#    L1NonReference.Rdata: results from our analysis

# Output:
#   
#    : ...

##############################################

#######################################
#                                     #
#    Organize eul1db data             #
#                                     #
#######################################

# Load necessary packages
library(rtracklayer)

# Read in eul1db data
Individuals <- read.delim("D:/HumanGenome/Data/Individuals_eul1db", skip = 5)
Samples     <- read.delim("D:/HumanGenome/Data/Samples_eul1db", skip = 5)
MRIP        <- read.delim("D:/HumanGenome/Data/MRIP_eul1db", skip = 5)
SRIP        <- read.delim("D:/HumanGenome/Data/SRIP_eul1db", skip = 5)

# Get sample ID for NA12878
idxNA12878 <- which(Samples$Sample_name == "NA12878")
ID1878     <- Samples$Individual_id[Samples$Sample_name == "NA12878"][1]
Samples[Samples$Sample_name == "NA12878",]

############
# Find NA12878 entries in MRIP
############

# Get entries with that sample ID in MRIP and replace strand'.' by '*'
idxMRIP_NA12878 <- grep(ID1878, MRIP$samples)
MRIP_NA12878    <- MRIP[idxMRIP_NA12878, c("chromosome", "X.mrip_accession.no",
                     "start", "stop", "strand", "referencel1hs", "studies", 
                     "samples")]
MRIP_NA12878$strand <- as.character(MRIP_NA12878$strand)
MRIP_NA12878$strand[MRIP_NA12878$strand == "."] <- "*"
grep("Beck2010", MRIP_NA12878$studies)

# Get genomic ranges of NA12878 L1HS from MRIP
IR_NA12878 <- IRanges(start = MRIP_NA12878$start, end = MRIP_NA12878$stop)
GR_NA12878 <- GRanges(seqnames = as.character(MRIP_NA12878$chromosome),
                      ranges = IR_NA12878, strand = MRIP_NA12878$strand)

# Match coordinates to hg38
GR_NA12878_Hg38 <- liftOver(GR_NA12878, 
   chain = import.chain("D:/L1polymORF/Data/hg19ToHg38.over.chain"))

############
# Find NA12878 entries in SRIP
############

# Get entries with that sample ID in SRIP and replace strand'.' by '*'
idxSRIP_NA12878 <- which(SRIP$sampleID  %in% idxNA12878)
SRIP_NA12878    <- SRIP[idxSRIP_NA12878, c("chromosome", "integrity", "mrip", "study_id",
   "g_start", "g_stop", "ref.start", "ref_stop","g_strand","referencel1hs")]
SRIP_NA12878[SRIP_NA12878$integrity == "full-length",]
SRIP_NA12878$g_strand <- as.character(SRIP_NA12878$g_strand)
SRIP_NA12878$g_strand[SRIP_NA12878$g_strand == "."] <- "*"

# Get genomic ranges of NA12878 L1HS
IR_NA12878_S <- IRanges(start = SRIP_NA12878$g_start, end = SRIP_NA12878$g_stop)
GR_NA12878_S <- GRanges(
  seqnames = paste("ch", as.character(SRIP_NA12878$chromosome), sep = ""),
  ranges = IR_NA12878_S, strand = SRIP_NA12878$g_strand)

# Match coordinates to hg38
GR_NA12878_S_Hg38 <- liftOver(GR_NA12878_S, 
   chain = import.chain("D:/L1polymORF/Data/hg19ToHg38.over.chain"))

############
# Explore similarities between MRIP_NA12878 and SRIP_NA12878
############

# 
Munique <- as.character(unique(SRIP_NA12878$mrip))
sum(Munique  %in% MRIP_NA12878$X.mrip_accession.no)
sum(!Munique %in% MRIP_NA12878$X.mrip_accession.no)


countOverlaps(GR_NA12878_S, GR_NA12878)

#######################################
#                                     #
#    Compare to L1NonReference        #
#                                     #
#######################################

# Load results
load("D:/L1polymORF/Data/L1NonReference.Rdata")

# Increase size in ranges from euL1db
GR_NA12878_Hg38_Large <- resize(GR_NA12878_Hg38, 2, fix="center")

# Find overlaps between islands and L1HS ranges
blnOverlapIslands_All <- overlapsAny(IslGRanges_reduced, L1GRanges)
blnOverlapBetwData    <- overlapsAny(IslGRanges_reduced, GR_NA12878_Hg38)
blnOverlapBetwData    <- overlapsAny(L1GRanges, GR_NA12878_Hg38_Large)
sum(blnOverlapBetwData)
IslGRanges_nonRef <- IslGRanges_reduced[!blnOverlapIslands_All]
IslGRanges_inRef  <- IslGRanges_reduced[blnOverlapIslands_All]

blnInMRIP <- overlapsAny(IslGRanges_nonRef, GR_NA12878_Hg38)
sum(blnInMRIP)

blnInMRIP <- overlapsAny(IslGRanges_inRef, GR_NA12878_Hg38)
sum(blnInMRIP)

Dist1 <- distanceToNearest(unlist(GR_NA12878_Hg38), L1GRanges)
Dist2 <- distanceToNearest(unlist(GR_NA12878_Hg38), IslGRanges_inRef)
hist(Dist1@elementMetadata@listData$distance, breaks = seq(0, 41000, 500))
sum(Dist1@elementMetadata@listData$distance < 10000)

# Look at scanned ranges
NrMapped <- sapply(ScannedL1Ranges, function(x) sum(!is.na(x[[1]]$pos)))
which.max(NrMapped)
FileNames[which.max(NrMapped)]

#######################################
#                                     #
#    Get L1 coordinates from Beck2010        #
#                                     #
#######################################

# Read in file with sequence
list.files('D:/L1polymORF/Data/')
BeckData <- read.csv('D:/L1polymORF/Data/Beck2010_mergedTable_withEmptySite.csv')
library(BSgenome.Hsapiens.UCSC.hg19)

idxNA12878_Beck <- which(BeckData$Coriell_ID == "NA12878")

CutOff <- 16
PMRangeList <- lapply(idxNA12878_Beck, function(i){
  Chrom <- paste("chr", BeckData$Chromosome[i], sep = "")
  Pattern <- as.character(BeckData$EmptySite[i])
#  Pattern <- substr(Pattern, 1, CutOff)
  PM <- matchPattern(Pattern, BSgenome.Hsapiens.UCSC.hg19[[Chrom]])
  ranges(PM)
})
idxMatch <- which(sapply(PMRangeList, length) > 0)
PMGRangeList <- lapply(idxMatch, function(i){
  Chrom <- paste("chr", BeckData$Chromosome[idxNA12878_Beck[i]],sep = "")
  GRanges(seqnames = Chrom, ranges = PMRangeList[[i]])
})
PMGRangeList <- GRangesList(PMGRangeList)
GRangesBeck2010 <- unlist(PMGRangeList)
overlapsAny(GRangesBeck2010, GR_NA12878) 

GRangesBeck2010_Hg38 <- liftOver(GRangesBeck2010, 
   chain = import.chain("D:/L1polymORF/Data/hg19ToHg38.over.chain"))

overlapsAny(GRangesBeck2010, IslGRanges_inRef) 

