##############################################
#
# General description:
#
#   The following script reads a bam for 12878 and explores peak
#   overlap 

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(pROC)
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(QuasR)

# Get all ranges of reads for per chromosome
Chroms       <- paste('chr', c(1:22, "X", "Y"), sep = "")

# Files and folders
BamFile           <- "D:/L1polymORF/Data/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
InFastQFolder     <- "D:/L1polymORF/Data/FastQ"
OutFastQFolder    <- "D:/L1polymORF/Data/FastqPerSuspectPeak/"
SampleFileFolder  <- "D:/L1polymORF/Data/"
L1TableFileName   <- "D:/L1polymORF/Data/L1_repeat_table.csv"

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName)

# Create names of subfamilies
repNChar    <- as.character(L1Table$repName)
SubFamiliesLookUp <- sapply(unique(repNChar), function(x){
  c(Name = x,
    SubFam = paste(strsplit(x, '[1-9]')[[1]][1:2], collapse = "1"))
})
NameMatch <- match(repNChar, SubFamiliesLookUp["Name",])
SubFamilies <- SubFamiliesLookUp["SubFam", NameMatch]

# Create GRanges objects with L1 Seqences
L1IRanges <- IRanges(start = L1Table$genoStart,
                     end = L1Table$genoEnd)
L1GRanges <- GRanges(seqnames = L1Table$genoName, ranges = L1IRanges,
                     strand = L1Table$strand)
GRanges_L1HSFL <- L1GRanges[width(L1GRanges) >= 6000 & SubFamilies == "L1HS"]
GRanges_L1PAFL <- L1GRanges[width(L1GRanges) >= 5900 & SubFamilies == "L1PA"]

cat("*******   Turning BAM files into GRanges ...   *******\n")

# Read coverage per chromosome
CoverList <- lapply(Chroms, function(Chrom){
  cat("Reading reads for chromosome", Chrom, "\n")
  ChromLength <- length(BSgenome.Hsapiens.UCSC.hg19[[Chrom]])
  R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))
  Reads   <- extractReads(bam.file = BamFile , region = R1)
  ReadCov <- coverage(Reads)
})

#######################################
#                                     #
#    Determine 'islands' and          #
#        overlap with L1              #
#                                     #
#######################################

# Determine separate islands with continuous read coverage and turn islands into
# genomic ranges
IslandList <- lapply(CoverList, function(x){
  Islands <- slice(x, lower = 1)
})
IslandGRanges <- lapply(1:length(IslandList), function(i){
  GRanges(seqnames = Chroms[i], 
          ranges = IslandList[[i]]@listData[[1]]@ranges,
          coverTotal = viewSums(IslandList[[i]])[[1]],
          coverMax   = viewMaxs(IslandList[[i]])[[1]])
})
IslandGRanges <- GRangesList(IslandGRanges)
IslandGRanges <- unlist(IslandGRanges)
RangeHist1    <- hist(width(IslandGRanges), breaks = seq(0, 25000, 100), 
                   plot = F)

# Find overlaps between islands and L1HS ranges
blnOverlapIslands_All  <- overlapsAny(IslandGRanges, L1GRanges)
blnOverlapIslands_L1HS <- overlapsAny(IslandGRanges, GRanges_L1HSFL)
blnOverlapIslands_L1PA <- overlapsAny(IslandGRanges, GRanges_L1PAFL)
blnOverlapL1           <- overlapsAny(GRanges_L1HSFL, IslandGRanges)

# Subset to get islands overlapping with L1HS
IslandGRangesSubset <- IslandGRanges[blnOverlapIslands_L1HS]
RangeHist2 <- hist(width(IslandGRangesSubset), breaks = seq(0, 25000, 100), 
                   plot = F)
plot(RangeHist1$mids, RangeHist1$density, type = "l", ylab = "Density", 
     xlab = "Width [bp]", ylim = c(0, 0.0005))
lines(RangeHist2$mids, RangeHist2$density, col = "red")

#######################################
#                                     #
#    Explore peak calling             #
#                                     #
#######################################

# Determine receiver-operating curves based on total and maximum coverage 
cat("*******   Calculating ROC curves ...   *******\n")
# rocTotal <- pROC::roc(response = blnOverlapIslands_L1HS,
#                 predictor = IslandGRanges@elementMetadata@listData$coverTotal)
# rocMax   <- pROC::roc(response = blnOverlapIslands_L1HS,
#                 predictor = IslandGRanges@elementMetadata@listData$coverMax)
# 
# # Plot ROCs
# Cols <- rainbow(2)
# plot(1 - rocTotal$specificities, rocTotal$sensitivities,
#      xlab = "1 - specificity", ylab = "Sensitivity", col = Cols[1])
# points(1 - rocMax$specificities, rocMax$sensitivities, col = Cols[2])
# legend("bottomright", legend = c("Total coverage", "Maximum coverage"),
#        col = Cols, lty = c(1,1))
# CreateDisplayPdf('D:/L1polymORF/Figures/PeakCallingROC.pdf')
# 
# rocMax$thresholds[rocMax$specificities > 0.9 & rocMax$sensitivities > 0.6]
hist(IslandGRanges@elementMetadata@listData$coverMax,
     breaks = 0:2500, xlim = c(0, 300), ylim  = c(0, 1000))

#######################################
#                                     #
#    Explore peak calling             #
#                                     #
#######################################

# Subset to get islands overlapping with L1HS
maxCover      <- IslandGRanges@elementMetadata@listData$coverMax
maxCoverRatio <- IslandGRanges@elementMetadata@listData$coverMax / 
  width(IslandGRanges)
MaxHist1 <- hist(maxCover[blnOverlapIslands_All], breaks = seq(0, 2000, 10),
                 plot = F)
MaxHist2 <- hist(maxCover[!blnOverlapIslands_All], breaks = seq(0, 3000, 10),
                 ylim = c(0, 100))
MaxHist2 <- hist(maxCover[!blnOverlapIslands_All], breaks = seq(0, 3000, 10),
                 plot = F)
plot(MaxHist1$mids, MaxHist1$counts, type = "s", ylab = "Density", 
     xlab = "Width [bp]", ylim = c(0, 1000), xlim = c(0, 1000))
lines(MaxHist2$mids, MaxHist2$counts, col = "red", type = "s")
MaxHist1$density[1:70]

sum(maxCover[!blnOverlapIslands_All] > 200)
sum(!blnOverlapIslands_All)
sum(blnOverlapIslands_All)

# RatioHist1 <- hist(maxCoverRatio[blnOverlapIslands_All], breaks = seq(0, 2, 0.01))
# RatioHist2 <- hist(maxCoverRatio[!blnOverlapIslands_All], breaks = seq(0, 5, 0.01),
#                  ylim = c(0, 10))
# RatioHist2 <- hist(maxCoverRatio[!blnOverlapIslands_All], breaks = seq(0, 3000, 10),
#                  plot = F)
# plot(RatioHist1$mids, RatioHist1$counts, type = "s", ylab = "Density", 
#      xlab = "Width [bp]", ylim = c(0, 1000), xlim = c(0, 1000))
# lines(RatioHist2$mids, RatioHist2$counts, col = "red", type = "s")

#######################################################
#                                                     #
#    Write fastq of suspected L1 not in reference     #
#                                                     #
#######################################################

# Get ranges of suspected L1s
idxSuspectL1Ranges <- which(maxCover > 100 & (!blnOverlapIslands_All))
SuspectL1Ranges    <- IslandGRanges[idxSuspectL1Ranges]
nearest(SuspectL1Ranges, L1GRanges)

# Create a vector of fastq file names
FastQnames <- paste(as.vector(seqnames(SuspectL1Ranges)), 
                    idxSuspectL1Ranges, sep = "_")
FastQPaths <- paste(OutFastQFolder, FastQnames, ".fastq", sep = "")
if (!dir.exists(OutFastQFolder)){
  dir.create(OutFastQFolder)
}
SampleFilePaths <- paste(SampleFileFolder, FastQnames, "Sample.txt", 
                         sep = "")

# Save sample file tables to be used by qAlign
SampleFileLine1 <-  c("FileName",	"SampleName")
for (i in 1:length(SampleFilePaths)){
  FileTable <- rbind(SampleFileLine1, c(FastQPaths[i], FastQnames[i]))
  write.table(FileTable, SampleFilePaths[i], sep = "\t", quote = F, 
              row.names = F, col.names = F)
}

# Write little fastq files per suspected 
WriteFastQPerRange(Ranges = SuspectL1Ranges, 
                   InBamfilePath  = BamFile,
                   InFastQfilePaths = list.files(InFastQFolder, full.names = T),
                   OutFilePaths = FastQPaths) 

#######################################
#                                     #
#     Map fastq file per range        #
#            to L1HS                  #
#                                     #
#######################################
i <- 1
# for (i in 1:length(SampleFilePaths)){
#   cat("Mapping peak", i, "of",  length(SampleFilePaths), "\n")
#   qAlign(sampleFile = SampleFilePaths[i],
#          genome = "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa", 
#          paired = "no")
# }

#######################################
#                                     #
#     Import reads mapped to L1       #
#                                     #
#######################################

# Get all bam files of reads mapped to L1HS
# FileNames <- list.files(OutFastQFolder, pattern = "chr[1-9]_[1-9]",
#                         full.names = T)
# FileNames <- FileNames[grep(".bam", FileNames)]
# FileNames <- FileNames[-grep(".bam.", FileNames)]
SamFileNames <- list.files(OutFastQFolder, pattern = "_aln.sam",
                        full.names = T)
fn <- SamFileNames[1]
for (fn in SamFileNames) {
  asBam(fn, destination = substr(fn, 1, nchar(fn) - 4))
}
FileNames <- list.files(OutFastQFolder, pattern = "_aln.bam",
                           full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]



# Loop through file names and read in bam files of reads mapped to L1
ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})

# Plot histogram of number mapped
hist(NrMapped2L1, breaks = 0:120)
sum(NrMapped2L1 > 0)

# Get aligned reads per peak
R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", ranges = IRanges(start = 1, end = 6000))
ReadsPerL1 <- lapply(FileNames[NrMapped2L1 > 0], function(x) {
  Reads <- extractReads(x, R1)
})

# Calculate a coverage matrix
CoverMat <- t(sapply(ReadsPerL1, function(x){
  Cov <- coverage(x)
  as.vector(Cov$L1HS_L1_Homo_sapiens)
}))

# Get means and 95% quantiles 
QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
idxFw <- 1:ncol(CoverMat)
idxRv <- ncol(CoverMat):1
plot(QuantileMat[2,], type = "n", ylim = c(0, 30), 
     ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)
CreateDisplayPdf('D:/L1polymORF/Figures/L1HSCoverNonReference.pdf')

# Determine range index from file name 
idxRange <- sapply(FileNames, function(x) as.numeric(strsplit(x, "_")[[1]][2]))
plot(maxCover[idxRange], NrMapped2L1, xlab = "Maximum coverage in range", 
     ylab = "Number reads mapped to L1HS")
CreateDisplayPdf('D:/L1polymORF/Figures/Range_Vs_L1HSCoverage.pdf')

