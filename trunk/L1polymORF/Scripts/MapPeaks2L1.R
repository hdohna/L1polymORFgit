##############################################
#
# General description:
#
#   The following script reads bam files of reads coming from chip-seq with
#   capture oligos containing L1HS sequences. It loops through the peaks,
#   get all reads from a peak and maps them on a L1HS sequence

# Input:
#
#    SRIP_eul1db: data on methods

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
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QuasR)

# Get all ranges of reads for per chromosome
Chromosomes <- paste("chr", 1:22, sep = "")

# Files and folders
DataFolder <- "D:/L1polymORF/Data/"
L1ReferenceBedFile <- "D:/L1polymORF/Data/hg38.fa_L1HS_6kb.bed"
PeakBedFile <- "D:/L1polymORF/Data/NA12878_L1_capt_IDTx_peaks.bed"
PeakBamFile <- "D:/L1polymORF/Data/NA12878.recal.sorted.bam"
FastQFile   <- "D:/L1polymORF/Data/NA12878.recal.sorted.bam"
FastQFolder <- "D:/L1polymORF/Data/L1HS cap-22668675/LHS0001trmIDXRCNMP3NMP3Q20430-25326778"
list.files("D:/L1polymORF/Data/L1HS cap-22668675/LHS0001trmIDXRCNMP3NMP3Q20430-25326778")

########################################
#                                      #
#        Define functions              #
#                                      #
########################################

# Function to write reads from a peak as fastq file and create sample file for
# the function qAlign
WriteFastqAndSample <- function(ReadList, 
                                Folder = DataFolder,
                                FilePrefix = "Peak"){
  # Create file name and path for fastq file
  FileNameFastq <- paste(FilePrefix, "fastq", sep = ".")
  FilePathFastq <- paste(Folder, FileNameFastq, sep = "")
  
  # Create lines of the fastq file and save them
  Seqs   <- as.character(ReadList$seq)
  Quals  <- as.character(ReadList$qual)
  SeqIDs <- paste(ReadList$rname, ReadList$pos, sep = ":")
  FastQLines1 <- paste("@", SeqIDs, sep = "")
  FastQLines3 <- paste("+", SeqIDs, sep = "")
  FastQLinesAll <- rep(NA, 4 * length(FastQLines1))
  FastQLinesAll[seq(1, length(FastQLinesAll) - 3, 4)] <- FastQLines1
  FastQLinesAll[seq(2, length(FastQLinesAll) - 2, 4)] <- Seqs
  FastQLinesAll[seq(3, length(FastQLinesAll) - 1, 4)] <- FastQLines3
  FastQLinesAll[seq(4, length(FastQLinesAll) - 0, 4)] <- Quals
  writeLines(FastQLinesAll, FilePathFastq)
  
  # Create lines of the sample file (for function qAlign) and save them
  FileTable <- rbind(c("FileName",	"SampleName"), 
                     c(FileNameFastq, FilePrefix))
  
  # Create file name and path for fastq file
  FileNameSample <- paste(FilePrefix, "Sample", sep = "")
  FilePathSample <- paste(Folder, FileNameSample, sep = "")
  write.table(FileTable, FilePathSample, sep = "\t", quote = F, 
              row.names = F, col.names = F)
  FilePathSample
}

#######################################
#                                     #
#     Read peaks and                  #
#   save fastq file per peak          #
#                                     #
#######################################

# Import BAM files as GRanges      
cat("*******   Importing BED files as GRanges ...   *******\n")
L1Ref        <- import.bed(con = L1ReferenceBedFile) 
NA12878Peaks <- import.bed(con = PeakBedFile) 

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
NrReads <- sapply(ScannedRanges, function(x) length(x$qname))

# Get index with minimum read threshold
idxBigPeak <- which(NrReads > 500)

# Create a fastq and sample file names
FilePrefixes <- paste(Chrom, "PeakPair", idxBigPeak, sep = "_")
FileNameMat  <- sapply(c("SampleFile", "1.fastq", "2.fastq", ""), function(x){
  paste(FilePrefixes, x, sep = "_")
})
FilePathMat <- matrix(paste(DataFolder,FileNameMat[,-4], sep = ""), ncol = 3)
      
# Save sample file tables to be used by qAlign
SampleFileLine1 <-  c("FileName1",	"FileName2", "SampleName")
for (i in 1:length(FQFileNames1)){
  FileTable <- rbind(SampleFileLine1, FileNameMat[i, c(2: 4)])
  write.table(FileTable, FilePathMat[i,1], sep = "\t", quote = F, 
              row.names = F, col.names = F)
}

# Open one connection per fastq file
FastQFileNames <- list.files(FastQFolder, full.names = T)
NrReadsPerIter <- 10^6
Stream1 <- FastqStreamer(FastQFileNames[1], NrReadsPerIter)
Stream2 <- FastqStreamer(FastQFileNames[2], NrReadsPerIter)
NrReadsRead <- 1

# Loop through fastq file and append to little peak-specific fastq files
while (NrReadsRead > 0){
  Reads1  <- yield(Stream1)
  Reads2  <- yield(Stream2)
  NrReadsRead <- length(Reads1)
  
  # Loop over peaks and create a fastq file per peak
  cat("Looping over peaks and writing to fastq files per peak\n")
  if (NrReadsRead > 0) {
    for (i in 1:length(idxBigPeak)){
      
      # Get names of reads in current peak, subset reads of current chunk and 
      # append to fastq files
      idxPeak <- idxBigPeak[i]
      ReadNames   <- ScannedRanges[[idxPeak]]$qname
      ReadSubset  <- Reads1@id %in% ReadNames
      Reads1Local <- Reads1[ReadSubset]
      Reads2Local <- Reads2[ReadSubset]
      writeFastq(Reads1Local, FilePathMat[i, 2], mode = "a", compress = F) 
      writeFastq(Reads2Local, FilePathMat[i, 3], mode = "a", compress = F) 
    }
  }
}
close(Stream1)
close(Stream2)

#######################################
#                                     #
#     Read fastq file per peak        #
#            to L1HS                  #
#                                     #
#######################################

for (i in 1:length(idxBigPeak)){
cat("Mapping peak", i, "of",  length(idxBigPeak), "\n")
  qAlign(sampleFile = FilePathMat[i, 1],
         genome = "D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa", 
         paired = "fr")
}

# for (i in 1:length(ScannedRanges)){
#   FilePrefix <- paste(Chrom, "Peak", i, sep = "_")
#   sampleFile <- WriteFastqAndSample(ScannedRange[[i]], FilePrefix = FilePrefix)
# cat("Mapping peak", i, "of",  length(ScannedRanges), "\n")
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
FileNames <- list.files(DataFolder, pattern = "chr1_Peak_[1-9]",
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

# Plot number mapped vs total number per 
NrReads_ordered <- NrReads[order(as.character(1:331))]
plot(NrReads_ordered, NrMapped2L1)
plot( NrReads_ordered, NrMapped2L1/NrReads_ordered)

scanBamHeader(FileNames)

R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", ranges = IRanges(start = 1, end = 6000))
ReadsPerL1 <- lapply(FileNames[NrMapped2L1 > 100], function(x) {
  Reads <- extractReads(x, R1)
})
ReadsPerL1[[3]]
i <- 3

# Calculate a coverage matrix
CoverMat <- t(sapply(ReadsPerL1, function(x){
  Cov <- coverage(x)
  as.vector(Cov$L1HS_L1_Homo_sapiens)
}))

Cov <- coverage(ReadsPerL1[[1]])
plot(colMeans(CoverMat), type = "s")
for (i in 1:length(ReadsPerL1)){
  Cov <- coverage(ReadsPerL1[[i]])
  lines(as.vector(Cov$L1HS_L1_Homo_sapiens), type = "s")
  
}
length(ReadsPerL1)
