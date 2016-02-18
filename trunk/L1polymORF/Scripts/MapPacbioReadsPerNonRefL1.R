##############################################
#
# General description:
#
#   The following script writes fastq files from pacbio reads of peaks
#   outside known L1 and aligns them to the reference L1

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1HSTableFileName: path to file that contains L1HS ranges in a table


# Output:
#   
#    ReadsPerNonRefL1.RData: ...

##############################################

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.R')

OutFastQFolder    <- "/home/hzudohna/L1polymORF/Data/PacbioFastqPerSuspectPeak/"
L1Consensus       <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
CoverSummaryPlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverNonReference_Pacbio.pdf'
CoverComparePlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverComparison_Pacbio.pdf'
OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference_Pacbio.Rdata'

# Suffices for alignment files created by BWA
SamSuffix <- "_aln.sam"
BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")

# BWA command (options can be added here)
BWAcommand <- '/home/txw/bwa/bwa-0.7.12/bwa mem'

# Load reads intersecting with suspected L1 ranges
load("/home/hzudohna/L1polymORF/Data/ReadsPerNonRefL1.RData")

#######################################################
#                                                     #
#    Write fastq of suspected L1 not in reference     #
#                                                     #
#######################################################

cat("*******   Writing little fastq files ...   *******\n")

# creat a vector of file prefixes
FilePrefix <- names(ScannedReads)
FilePrefix <- gsub(":", "_", FilePrefix)

FilePaths <- t(sapply (1:length(ScannedReads), function(i){
  cat("Writing reads for peak", i, "of", length(ScannedReads), "\n")
  Reads <- ScannedReads[[i]]
  WriteFastqAndSample(Reads, OutFastQFolder, FilePrefix[i])
}))

FastQPaths <- FilePaths[,"FilePathFastq"]

#######################################
#                                     #
#     Map fastq file per range        #
#            to L1HS                  #
#                                     #
#######################################

cat("*******   Mapping little fastqs to L1 ...   *******\n")

# Create index file
CmdIndex <- paste('bwa index', L1Consensus)
system(CmdIndex)

# Run BWA for each little fastq file  
OutFiles <- paste(substr(FastQPaths, 1, nchar(FastQPaths) - 6), SamSuffix, sep = "")
CmdLines <- paste(BWAcommand,  L1Consensus, FastQPaths)
CmdLines <- paste(CmdLines, OutFiles, sep = " > ")
for (CmdL in CmdLines) system(CmdL)

#######################################
#                                     #
#     Import reads mapped to L1       #
#                                     #
#######################################

cat("*******  Reading and analyzing mapped reads ...   *******\n")

# Get all names of sam files created by BWA
SamFileNames <- list.files(OutFastQFolder, pattern = SamSuffix,
                           full.names = T)

# Turn sam files into bam files
for (fn in SamFileNames) {
  asBam(fn, destination = substr(fn, 1, nchar(fn) - 4))
}

# get names of newly created bam files
FileNames <- list.files(OutFastQFolder, pattern = BamSuffix,
                        full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Loop through file names and read in bam files of reads mapped to L1
ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})

# Get aligned reads per peak
R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
              ranges = IRanges(start = 1, end = 6000))
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
pdf(file = CoverSummaryPlot)
plot(QuantileMat[2,], type = "n", ylim = c(0, max(QuantileMat)), 
     ylab = 'Coverage', xlab = "Genomic position")
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,], lwd = 1.2)
dev.off()

# Determine range index from file name 
pdf(file = CoverComparePlot)
idxRange <- sapply(FileNames, function(x) as.numeric(strsplit(x, "_")[[1]][2]))
plot(maxCover[idxRange], NrMapped2L1, xlab = "Maximum coverage in range", 
     ylab = "Number reads mapped to L1HS")
dev.off()

# Save results
cat("*******  Saving results ...   *******\n")
save(list = c("FileNames", "ScannedL1Ranges", "ReadsPerL1", 
              "CoverMat", "QuantileMat", "idxRange"), file = OutResults)
