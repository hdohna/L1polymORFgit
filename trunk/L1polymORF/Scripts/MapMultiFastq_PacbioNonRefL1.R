##############################################
#
# General description:
#
#   The following script aligns all fastq files in a folder
#   to the same reference using bwa

# Input:
#
#     FastQFolder: folder that contains the fastq files to be mapped
#     Homo_sapiens_L1_consensus.fa: L1 consensus sequence

# Output:
#   
#    Little sam file for each fastq file

##############################################

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load necessary libraries
library(GenomicRanges)
library(Rsamtools)
library(csaw)

# Specify file paths 
OutFastQFolder    <- "/home/hzudohna/L1polymORF/Data/PacbioFastqPerSuspectPeak/"
L1Consensus       <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
CoverSummaryPlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverNonReference_Pacbio.pdf'
CoverComparePlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverExamples_Pacbio.pdf'
OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference_Pacbio.Rdata'

# Suffices for alignment files created by BWA
SamSuffix <- "_aln.sam"
BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")

#######################################
#                                     #
#     Map fastq file per range        #
#            to L1HS                  #
#                                     #
#######################################

# # Map all fastq files in FastQFolder to L1 consensus 
# FilePaths <- MapMultiFastq(FastQFolder = OutFastQFolder,
#    AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0',
#    Reference = L1Consensus)
                          
#######################################
#                                     #
#     Import reads mapped to L1       #
#                                     #
#######################################

cat("*******  Reading and analyzing mapped reads ...   *******\n")

# Get all names of sam files created by BWA
#SamFileNames <- FilePaths$SamPath

# Get all names of sam files created by BWA
SamFileNames <- list.files(OutFastQFolder, pattern = SamSuffix,
                           full.names = T)

# Turn sam files into bam files
for (fn in SamFileNames) {
  asBam(fn, destination = substr(fn, 1, nchar(fn) - 4), overwrite = T)
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

# # Determine range index from file name 
pdf(file = CoverComparePlot)
plot(CoverMat[1,], type = "s", xlab = "Position on L1",
     ylab = "Coverage", ylim = c(0, 100))
Cols <- rainbow(nrow(CoverMat))
for (i in 1:nrow(CoverMat)){
  lines(CoverMat[i,], type = "s", col = Cols[i])
  
}
dev.off()

# Save results
cat("*******  Saving results ...   *******\n")
save(list = c("FileNames", "ScannedL1Ranges", "ReadsPerL1", 
              "CoverMat", "QuantileMat"), file = OutResults)
