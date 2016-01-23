##############################################
#
# General description:
#
#   The following script reads bam files of reads coming from chip-seq with
#   capture oligos containing L1HS sequences. It loops through the peaks,
#   get all reads from a peak and maps them on a L1HS sequence

# Input:
#
#     DataFolder: path to folder where all output data should be saved
#     L1ReferenceBedFile: path to file that contains genomic ranges of L1
#     PeakBedFile: path to file that contains genomic ranges of mapped peaks
#     PeakBamFile: path to file that contains mapped reads
#     FastQFolder: path to folder that contains fastq files
#     L1HSTableFileName: path to file that contains L1HS ranges in a table


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
library(ShortRead)
library(csaw)
library(chipseq)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QuasR)

# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")

# Define a minimum number of reads per peak for it to be mapped to L1
MinReadPerPeak <- 500

# Files and folders
DataFolder <- "D:/L1polymORF/Data/"
L1ReferenceBedFile <- "D:/L1polymORF/Data/hg38.fa_L1HS_6kb.bed"
PeakBedFile <- "D:/L1polymORF/Data/NA12878_L1_capt_IDTx_peaks.bed"
PeakBamFile <- "D:/L1polymORF/Data/NA12878.recal.sorted.bam"
#FastQFolder <- "D:/L1polymORF/Data/L1HS cap-22668675/LHS0001trmIDXRCNMP3NMP3Q20430-25326778"
FastQFolder <- "D:/L1polymORF/Data/FastQ"
L1HSTableFileName <- "D:/L1polymORF/Data/L1HS_repeat_table.csv"

# Load functions

#######################################
#                                     #
#     Read peaks and                  #
#   save fastq file per peak          #
#                                     #
#######################################

# Import BED files as GRanges      
cat("*******   Importing BED files as GRanges ...   *******\n")
L1Ref        <- import.bed(con = L1ReferenceBedFile) 
NA12878Peaks <- import.bed(con = PeakBedFile)
ranges(NA12878Peaks)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRangesBED  <- IRanges(start = L1Ref@elementMetadata$thick@width,
                     end = L1Ref@elementMetadata$thick@width +
                       L1Ref@elementMetadata$thick@start)
L1GRangesBED  <- GRanges(seqnames = seqnames(L1Ref), ranges = L1IRangesBED ,
                     strand = strand(L1Ref))

# Read in table with L1 ranges
L1HSTable <- read.csv(L1HSTableFileName)

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = L1HSTable$genoStart,
                     end = L1HSTable$genoEnd)
L1GRanges <- GRanges(seqnames = L1HSTable$genoName, ranges = L1IRanges,
                     strand = L1HSTable$strand)
L1GRanges <- L1GRanges[width(L1GRanges) >= 6000]

# Get all peaks that overlap with L1 ranges on the reference genome
OverlapsWithL1HS <- findOverlaps(NA12878Peaks, L1GRanges)

# Get all L1 sequences  
L1Seq      <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GRanges, as.character = T)
L1SeqNames <- paste(as.vector(seqnames(L1GRanges)), start(L1GRanges), end(L1GRanges),
                    strand(L1GRanges), sep = "_")
# Loop over chromosomes      
cat("*******   Looping over chromosomes ...   *******\n")
Chrom <- "chr1"
cat("Analyzing chromosome", Chrom, "...\n")

cat("Reading reads\n")
RangesToScan <- RangesList(ranges(NA12878Peaks[seqnames(NA12878Peaks) == Chrom]))
names(RangesToScan) <- Chrom
param <- ScanBamParam(which = RangesToScan, what = scanBamWhat())
ScannedRanges <- scanBam(file = PeakBamFile, 
                        param = param)

# Get index with minimum read threshold
NrReads    <- sapply(ScannedRanges, function(x) length(x$qname))
idxBigPeak <- which(NrReads > MinReadPerPeak)

# Create a vector of all read names and peak indices
AllReadNames <- unlist(lapply(idxBigPeak, function(x){
  ScannedRanges[[x]]$qname
})) 
AllReadIndices <- unlist(lapply(1:length(idxBigPeak), function(x){
  rep(x, length(ScannedRanges[[idxBigPeak[x]]]$qname))
})) 
ReadNameDupl   <- duplicated(AllReadNames)
AllReadNames   <- AllReadNames[!ReadNameDupl]
AllReadIndices <- AllReadIndices[!ReadNameDupl]

# Create a fastq and sample file names
FilePrefixes <- paste(Chrom, "PeakFQ", idxBigPeak, sep = "_")
FileNameMat  <- sapply(c("SampleFile", ".fastq", ""), function(x){
  paste(FilePrefixes, x, sep = "_")
})
FileNameMat <- matrix(FileNameMat, ncol = 3)
FilePathMat <- matrix(paste(DataFolder,FileNameMat[,-3], sep = ""), ncol = 2)
      
# Save sample file tables to be used by qAlign
SampleFileLine1 <-  c("FileName",	"SampleName")
for (i in 1:length(FilePrefixes)){
  FileTable <- rbind(SampleFileLine1, FileNameMat[i, c(2:3)])
  write.table(FileTable, FilePathMat[i,1], sep = "\t", quote = F, 
              row.names = F, col.names = F)
}

# # Create a fastq and sample file names
# FilePrefixes <- paste(Chrom, "PeakPair", idxBigPeak, sep = "_")
# FileNameMat  <- sapply(c("SampleFile", "1.fastq", "2.fastq", ""), function(x){
#   paste(FilePrefixes, x, sep = "_")
# })
# FilePathMat <- matrix(paste(DataFolder,FileNameMat[,-4], sep = ""), ncol = 3)
#       
# # Save sample file tables to be used by qAlign
# SampleFileLine1 <-  c("FileName1",	"FileName2", "SampleName")
# for (i in 1:length(FilePrefixes)){
#   FileTable <- rbind(SampleFileLine1, FileNameMat[i, c(2:4)])
#   write.table(FileTable, FilePathMat[i,1], sep = "\t", quote = F, 
#               row.names = F, col.names = F)
# }

# Open one connection per fastq file
FastQFileNames <- list.files(FastQFolder, full.names = T)
NrReadsPerIter <- 10^6
Stream1 <- FastqStreamer(FastQFileNames[1], NrReadsPerIter)
Stream2 <- FastqStreamer(FastQFileNames[2], NrReadsPerIter)
NrReadsRead <- 1
ReadCounter <- 0
WriteCounter <- 0

# Loop through fastq file and append to little peak-specific fastq files
while (NrReadsRead > 0){
  Reads1  <- yield(Stream1)
  Reads2  <- yield(Stream2)
  NrReadsRead  <- length(Reads1)
  ReadCounter  <- ReadCounter + NrReadsRead
  ReadIDChar   <- as.character(Reads1@id)
  ReadIDsuffix <- substr(ReadIDChar, nchar(ReadIDChar) - 7, nchar(ReadIDChar))
  unique(ReadIDsuffix)
  ReadIDs      <- substr(ReadIDChar, 1, nchar(ReadIDChar) - 8)
  ReadSubset   <- ReadIDs %in% AllReadNames
  ReadIDs      <- ReadIDs[ReadSubset]
  Reads1       <- Reads1[ReadSubset]
  Reads2       <- Reads2[ReadSubset]
  WriteCounter <- WriteCounter + length(ReadIDs)
  cat("Total of", ReadCounter, "reads read \n")
  # Loop over peaks and create a fastq file per peak
  cat("Looping over peaks and writing to fastq files per peak\n")
  if (sum(ReadSubset) > 0) {
    Indices <- unique(AllReadIndices[AllReadNames %in% ReadIDs])
    for (i in Indices){
      
      # Get names of reads in current peak, subset reads of current chunk and 
      # append to fastq files
      idxPeak <- idxBigPeak[i]
      ReadNames   <- ScannedRanges[[idxPeak]]$qname
      ReadSubset  <- ReadIDs %in% ReadNames
      Reads1Local <- Reads1[ReadSubset]
      Reads2Local <- Reads2[ReadSubset]
      writeFastq(Reads1Local, FilePathMat[i, 2], mode = "a", compress = F) 
      writeFastq(Reads2Local, FilePathMat[i, 2], mode = "a", compress = F) 
    }
  }
  cat("Total of", WriteCounter, "written out \n")
}
close(Stream1)
close(Stream2)

#######################################
#                                     #
#     Map fastq file per peak         #
#            to L1HS                  #
#                                     #
#######################################

for (i in 1:length(idxBigPeak)){
cat("Mapping peak", i, "of",  length(idxBigPeak), "\n")
  qAlign(sampleFile = FilePathMat[i, 1],
         genome = "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa", 
         paired = "no")
}

# for (i in 1:length(idxBigPeak)){
#   j <- idxBigPeak[i]   
#   FilePrefix <- paste(Chrom, "Peak", j, sep = "_")
#   sampleFile <- WriteFastqAndSample(ScannedRanges[[j]], FilePrefix = FilePrefix)
# cat("Mapping peak", i, "of",  length(idxBigPeak), "\n")
#   qAlign(sampleFile = sampleFile,
#          genome = "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa", 
#          paired = "no")
# }

#######################################
#                                     #
#     Import reads mapped to L1       #
#                                     #
#######################################

# Get all bam files of reads mapped to L1HS
FileNames <- list.files(DataFolder, pattern = "chr1_PeakFQ_[1-9]",
                        full.names = T)
FileNames <- FileNames[grep(".bam", FileNames)]
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Loop through file names and read in bam files of reads mapped to L1
ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})

# Plot histogram of number mapped
hist(NrMapped2L1, breaks = 0:1600)
hist(NrMapped2L1, breaks = 0:1600, xlim = c(0, 100))
max(NrMapped2L1)

# Plot number mapped vs total number per peak
NrReads_ordered <- NrReads[order(as.character(1:331))]
plot(NrReads_ordered, NrMapped2L1)
plot( NrReads_ordered, NrMapped2L1/NrReads_ordered)

# Get aligned reads per peak
R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", ranges = IRanges(start = 1, end = 6000))
ReadsPerL1 <- lapply(FileNames[NrMapped2L1 > 10], function(x) {
  Reads <- extractReads(x, R1)
})

# Calculate a coverage matrix
CoverMat <- t(sapply(ReadsPerL1, function(x){
  Cov <- coverage(x)
  as.vector(Cov$L1HS_L1_Homo_sapiens)
}))

plot(colMeans(CoverMat), type = "s")
for (i in 1:length(ReadsPerL1)){
  Cov <- coverage(ReadsPerL1[[i]])
  lines(as.vector(Cov$L1HS_L1_Homo_sapiens), type = "s")
  
}

# Loop through file names get consensus sequence per file
NullSeq    <- rep('-', 6050)
ConsensSeq <- t(sapply(FileNames[NrMapped2L1 >= 100], function(x) {
  PU <- pileup(x)
  ConsensFromPileup(PU, NullSeq)
}))

# Get Peak IDs from consensus sequences
PeakIDsFromConsens <- as.numeric(sapply(rownames(ConsensSeq), 
                             function(x) strsplit(x, "_")[[1]][3]))
HitMatch <- match(PeakIDsFromConsens, OverlapsWithL1HS@queryHits)
L1Seq[OverlapsWithL1HS@subjectHits[HitMatch]]
which(!is.na(HitMatch))

HalfRange <- 10
Range     <- -HalfRange:HalfRange
sapply(which(!is.na(HitMatch)), function(i){
  ConSeq <- ConsensSeq[i,]
  SeqPos <- which(ConSeq != "-")
  L1DChar <-L1Seq[OverlapsWithL1HS@subjectHits[HitMatch[i]]]
#  L1DChar <-L1Seq[1]
  sum(s2c(L1DChar)[SeqPos[1:(length(SeqPos) - 2)]] != ConSeq[SeqPos[3:length(SeqPos)]])
  s2c(L1DChar)[SeqPos[1:10]]
  ConSeq[SeqPos[3:12]]
  NrFit <- sapply(Range, function(x) {
    sum(s2c(L1DChar)[SeqPos[(HalfRange + 1):(length(SeqPos) - HalfRange)]] == 
        ConSeq[SeqPos[(HalfRange + 1 + x):(length(SeqPos) - HalfRange + x)]])})
  OffSet <- Range[which.max(NrFit)]
  Start1 <- max(1, 1 + OffSet)
  End1   <- min(length(SeqPos), length(SeqPos) + OffSet)
  Start2 <- max(1, 1 - OffSet)
  End2   <- min(length(SeqPos), length(SeqPos) - OffSet)
  sum(ConSeq[SeqPos[Start1:End1]] != s2c(L1DChar)[SeqPos[Start2:End2]]) / length(SeqPos)
})
