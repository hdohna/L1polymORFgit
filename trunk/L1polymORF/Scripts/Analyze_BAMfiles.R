##############################################
#
# General description:
#
#   The following script reads bam files of reads coming from chip-seq with
#   capture oligos containing L1HS sequences

# Input:
#
#    SRIP_eul1db: data on methods

# Output:
#   
#    : ...

##############################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(VennDiagram)


# Read data on retrotranspsoson polymorphisms and put in in GRanges object
cat("*******   Reading SRIP data ...   *******\n")
SRIP <- read.delim("D:/L1polymORF/Data/SRIP_eul1db", skip = 5)
SRIPGRanges <- GRanges(seqnames = paste("chr", SRIP$chromosome, sep = ""),
   ranges = IRanges(start = pmin(SRIP$g_start,SRIP$g_stop),
                    end = pmax(SRIP$g_start,SRIP$g_stop)))
DuplStart <- duplicated(SRIP$g_start)
DuplEnd   <- duplicated(SRIP$g_stop)
AllowedStudies <- SRIP$study_id %in% c("Ewing2010", "Solyom2012")
AllowedStudies <- SRIP$study_id %in% unique(SRIP$study_id)
Subset <- (!DuplStart) & (!DuplEnd) & AllowedStudies
SRIP_NoDupl <- SRIP[Subset,] 
SRIPGRanges_NoDupl <- SRIPGRanges[Subset]

# Look at general patterns in the data
table(SRIP_NoDupl$study_id, SRIP_NoDupl$integrity)
table(SRIP_NoDupl$study_id, SRIP_NoDupl$lineage)
unique(SRIP_NoDupl$Sample_name)

# Read sample file
Samples     <- read.delim("D:/L1polymORF/Data/Samples_eul1db", skip = 5)
Samples[Samples$Sample_name %in% c("NA12878", "NA12892"), ]

#QA <- qa("D:/L1polymORF/Data/NA12878-L1HS_S1_L001_R1_001_13ac56fa66fa.bam",
#                             type ="BAM")
# 
#Al1 <- readAligned("D:/L1polymORF/Data/NA12878-L1HS_S1_L001_R1_001_13ac56fa66fa.bam",
#                  type = "BAM")


# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")
MinCoverage <- 5
IslRange    <- 10000
CalcNew_IslandPerCoverage <- T

#######################################
#                                     #
#    Turn BAM files into GRanges      #
#                                     #
#######################################

cat("*******   Turning BAM files into GRanges ...   *******\n")

# Function to get read coverage per chromosome
ReadCoverPerChrom <- function(InputFile) {
  lapply(Chromosomes, function(Chrom){
     ChromLength <- length(BSgenome.Hsapiens.UCSC.hg19[[Chrom]])
     R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))
     Reads <- extractReads(bam.file = InputFile, region = R1)
     ReadCov <- coverage(Reads)
   })
}

# Get read cover per chromosome for both alignments
if (file.exists("D:/L1polymORF/Data/GRanges_NA12878_NA12892.RData")){
  load(file = "D:/L1polymORF/Data/GRanges_NA12878_NA12892.RData")
} else {
  ReadCoverPerChrom1 <- 
    ReadCoverPerChrom("D:/L1polymORF/Data/NA12878-L1HS_S1_L001_R1_001_13ac56fa66fa.bam")
  ReadCoverPerChrom2 <- 
    ReadCoverPerChrom("D:/L1polymORF/Data/NA12892-L1HS_S2_L001_R1_001_13ac27cb1b4d.bam")
  save(file = "D:/L1polymORF/Data/GRanges_NA12878_NA12892.RData", 
       list = c("ReadCoverPerChrom1", "ReadCoverPerChrom2"))
}

#######################################
#                                     #
#      Get islands from GRanges       #
#                                     #
#######################################

cat("*******   Getting islands from GRanges ...   *******\n")

# Get overlaps with known SRIP per chromosome
IslandPerChrom <- function(ReadCovers, MinCoverage, IslRange){
  lapply(ReadCovers, function(ReadCov){
    Islands <- slice(ReadCov, lower = MinCoverage)
    IslandRanges <- resize(ranges(Islands), width = IslRange, fix = "center") 
    IslandRanges <- intersect(IslandRanges, IslandRanges) 
  })
}

# Determine the islands for different minimum coverage values
if (file.exists("D:/L1polymORF/Data/IslandPerCoverage.RData") & 
      (!CalcNew_IslandPerCoverage)){
  load(file = "D:/L1polymORF/Data/IslandPerCoverage.RData")
} else {
  MinCovVals <- 1:20
  IslandPerChrom1 <- lapply(MinCovVals, function(x) {
    IslandPerChrom(ReadCoverPerChrom1, x, IslRange)
  })
  IslandPerChrom2 <- lapply(MinCovVals, function(x) {
    IslandPerChrom(ReadCoverPerChrom2, x, IslRange)
  })
  save(file = "D:/L1polymORF/Data/IslandPerCoverage10Kbp.RData", 
       list = c("MinCovVals", "IslandPerChrom1", "IslandPerChrom2"))
}

#######################################
#                                     #
#        Analyzing islands            #
#                                     #
#######################################

cat("*******   Analyzing islands ...   *******\n")

cat("Count intersecting L1s ...\n")

# Determine per minimum coverage value the number of islands and the number of
# islands overlapping with SRIG
MeasuresPerMinCover <- function(IslandPerChrom) {
  sapply(IslandPerChrom, function(x){
     MeasuresPerChrom <- sapply(x, function(y) {
       c(Tot = length(y[[1]]), NrIntersect = sum(y %over% SRIPGRanges))
     })
     rowSums(MeasuresPerChrom)
   })
}
MeasuresPerMinCover1 <- MeasuresPerMinCover(IslandPerChrom1)
MeasuresPerMinCover2 <- MeasuresPerMinCover(IslandPerChrom2)

# Intersection between the two alignments
Intersect12 <- sapply(1:length(IslandPerChrom1), function(x){
  sum(sapply(1:length(IslandPerChrom1[[x]]), function(y){
    sum(IslandPerChrom1[[x]][[y]] %over% IslandPerChrom2[[x]][[y]])
  }))
})

# Function to get counts of interesecting L1s from the database separated
# by a factor (database column)
CountsPerFactor <- function(VariableName, IslandPerChrom){
  Vals   <- unique(SRIP[,VariableName])
  Counts <- rep(0, length(Vals))
  names(Counts) <- Vals
  CountsPerCover <- sapply(IslandPerChrom, function(x){
    CountsPerChrom <- sapply(x, function(y) {
      Subset <- SRIPGRanges_NoDupl %over% y
      NewCounts <- table(SRIP_NoDupl[Subset,VariableName])
      Counts[names(NewCounts)] <- NewCounts
      Counts
    })
    rowSums(CountsPerChrom)
  }) 
  SRIPCounts <- table(SRIP_NoDupl[,VariableName])
  SRIPCounts <- SRIPCounts[match(rownames(CountsPerCover), names(SRIPCounts))]
  cbind(CountsPerCover, SRIPCounts) 
}

# Determine the SRIG integrity per minimum coverage value  
# overlapping with the islands
cat("Count intersecting L1s per integrity level ...\n")
SRIGIntegrityPerMinCover1 <- CountsPerFactor("integrity", IslandPerChrom1)
SRIGIntegrityPerMinCover2 <- CountsPerFactor("integrity", IslandPerChrom2)

SRIPCounts <- SRIGIntegrityPerMinCover1[,"SRIPCounts"]
OverlapCounts1 <- SRIGIntegrityPerMinCover1[,1]
OverlapCounts2 <- SRIGIntegrityPerMinCover2[,1]

ConTable1 <- cbind(SRIPCounts - OverlapCounts1, OverlapCounts1)
colnames(ConTable1) <- c("NotDetected", "Detected")
chisq.test(ConTable1)
ConTable2 <- cbind(SRIPCounts - OverlapCounts2, OverlapCounts2)
colnames(ConTable2) <- c("NotDetected", "Detected")
chisq.test(ConTable2)

# Determine the SRIG subgroup per minimum coverage value  
# overlapping with the islands
# SRIGsub_groupPerMinCover1 <- CountsPerFactor("sub_group", IslandPerChrom1)
# SRIGsub_groupPerMinCover2 <- CountsPerFactor("sub_group",IslandPerChrom2)

# Determine the SRIG method per minimum coverage value  
# overlapping with the islands
cat("Count intersecting L1s per method ...\n")
SRIGMethodPerMinCover1 <- CountsPerFactor("method_id", IslandPerChrom1)
SRIGMethodPerMinCover2 <- CountsPerFactor("method_id", IslandPerChrom2)

# Determine the SRIG method per minimum coverage value  
# overlapping with the islands
cat("Count intersecting L1s per study ...\n")
SRIGstudy_idPerMinCover1 <- CountsPerFactor("study_id", IslandPerChrom1)
SRIGstudy_idPerMinCover2 <- CountsPerFactor("study_id", IslandPerChrom2)

#######################################
#                                     #
#         Plot results                #
#                                     #
#######################################

# plot the number of overlapping SRIPs, the total number of peaks
# and the proportion of peaks that intersect
par(mfrow = c(2, 2))
Y <- "NrIntersect"
plot(MinCovVals, MeasuresPerMinCover1[Y,], type = "l", 
     xlab = "Minimum coverage for peak calling",
     ylab = "Nr intersect with known SRIP", col = "red")
lines(MinCovVals, MeasuresPerMinCover2[Y,], col = "blue")
legend("topright", legend = c("NA12878", "NA12892"),
       col = c("red", "blue"), lty = rep(1, 2), cex = 0.5)

Y <- "Tot"
plot(MinCovVals, MeasuresPerMinCover1[Y,], type = "l", 
     xlab = "Minimum coverage for peak calling",
     ylab = "Total nr called", col = "red")
lines(MinCovVals, MeasuresPerMinCover2[Y,], col = "blue")
legend("topright", legend = c("NA12878", "NA12892"),
       col = c("red", "blue"), lty = rep(1, 2), cex = 0.5)

Y1 <- "NrIntersect"
Y2 <- "Tot"
plot(MinCovVals, MeasuresPerMinCover1[Y1,] / MeasuresPerMinCover1[Y2,], type = "l", 
     xlab = "Minimum coverage for peak calling",
     ylab = "Proportion intersect ", col = "red", ylim = c(0, 1))
lines(MinCovVals, MeasuresPerMinCover2[Y1,] / MeasuresPerMinCover2[Y2,], col = "blue")
lines(MinCovVals, Intersect12 / MeasuresPerMinCover2[Y2,], col = "green")
legend("topleft", legend = c("NA12878 & known SRIP", "NA12892 & known SRIP", "NA12878 & NA12892"),
       col = c("red", "blue", "green"), lty = rep(1, 3), cex = 0.5)

plot(MinCovVals,Intersect12, type = "l", 
     xlab = "Minimum coverage for peak calling",
     ylab = "Nr intersect betw. NA12878 & NA12892")

dev.copy2pdf(file = "D:/L1polymORF/Figures/PeakIntersection.pdf")

# Detect the proportion of L1 from the database that are detected by our method
# broken down by L1 integrity
par(mfrow = c(1, 1))
barplot(rbind(OverlapCounts1 / SRIPCounts, OverlapCounts2 / SRIPCounts),
        beside = T, legend.text = c("NA12878", "NA12892"), 
        names.arg = c("unknown", "5' truncated", "int. fragment", "full-length", "3' truncated"),
        args.legend = list(x = "topleft"),
        ylab = "Proportion detected", ylim = c(0, 0.07))
text(seq(2, 14, 3) - 0.5, OverlapCounts1 / SRIPCounts + 0.003, OverlapCounts1)
text(seq(2, 14, 3) + 0.5, OverlapCounts2 / SRIPCounts + 0.003, OverlapCounts2)
PropDetect1 <- sum(OverlapCounts1) / sum(SRIPCounts)
PropDetect2 <- sum(OverlapCounts2) / sum(SRIPCounts)
lines(c(0, 100), c(PropDetect1, PropDetect1), lty =2)
dev.copy2pdf(file = "D:/L1polymORF/Figures/DetectionByIntegrity10Kbp.pdf")

# Detect the proportion of L1 from the database that are detected by our method
# broken down by L1 integrity
par(mfrow = c(1, 1))
S1 <- rowSums(SRIGstudy_idPerMinCover1[,1:20]) > 0 & 
  rowSums(SRIGstudy_idPerMinCover2[,1:20]) > 0
PDetect1 <- SRIGstudy_idPerMinCover1[S1,1] /SRIGstudy_idPerMinCover1[S1,"SRIPCounts"]
PDetect2 <- SRIGstudy_idPerMinCover2[S1,1] /SRIGstudy_idPerMinCover2[S1,"SRIPCounts"]
barplot(rbind(PDetect1, PDetect2),
        beside = T, legend.text = c("NA12878", "NA12892"), 
        args.legend = list(x = "topleft"),
        ylab = "Proportion detected", ylim = c(0, 0.1),
        cex.names = 0.5)
text(seq(2, 23, 3) - 0.5, PDetect1 + 0.003, SRIGstudy_idPerMinCover1[S1,1])
text(seq(2, 23, 3) + 0.5, PDetect2 + 0.003, SRIGstudy_idPerMinCover2[S1,1])
dev.copy2pdf(file = "D:/L1polymORF/Figures/DetectionByStudy10Kbp.pdf")


