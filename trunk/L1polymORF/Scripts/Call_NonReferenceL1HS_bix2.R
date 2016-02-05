##############################################
#
# General description:
#
#   The following script reads a bam file, finds peaks that do not overlap
#   with reference L1s 

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
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get all ranges of reads for per chromosome
Chroms       <- paste('chr', c(1:22, "X", "Y"), sep = "")

# Files and folders
BamFile           <- "/home/hzudohna/BoData/NA12878-L1HS_S1_L001.dedup.unique.sorted.bam"
InFastQFolder     <- "/home/hzudohna/L1polymORF/Data/FastQ"
OutFastQFolder    <- "/home/hzudohna/L1polymORF/Data/FastqPerSuspectPeak/"
SampleFileFolder  <- "/home/hzudohna/L1polymORF/Data/"
L1TableFileName   <- "/home/hzudohna/L1polymORF/Data/L1_repeat_table.csv"
L1Consensus       <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
CoverSummaryPlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverNonReference.pdf'
CoverComparePlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverNonReference.pdf'
OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata'

# Suffices for alignment files created by BWA
SamSuffix <- "_aln.sam"
BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")

# BWA command (options can be added here)
BWAcommand <- 'bwa mem'

# Peak calling parameters
MinMaxCover <- 100    # minimum maximum coverage to be called a peak 
MinDist2L1  <- 3*10^4 # minimum distance to L1 to be called a peak 

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

cat("*******   Turning BAM files into GRanges ...   *******\n")

# Read coverage per chromosome
CoverList <- lapply(Chroms, function(Chrom){
  cat("Reading reads for chromosome", Chrom, "\n")
  ChromLength <- length(BSgenome.Hsapiens.UCSC.hg38[[Chrom]])
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

# Determine separate islands with continuous read coverage and turn islands 
# into genomic ranges
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

# Find overlaps between islands and L1HS ranges
blnOverlapIslands_All <- overlapsAny(IslandGRanges, L1GRanges)

#######################################################
#                                                     #
#    Write fastq of suspected L1 not in reference     #
#                                                     #
#######################################################

# Get ranges of suspected L1s
maxCover      <- IslandGRanges@elementMetadata@listData$coverMax
idxSuspectL1Ranges <- which(maxCover > MinMaxCover & (!blnOverlapIslands_All))
SuspectL1Ranges    <- IslandGRanges[idxSuspectL1Ranges]

# Remove ranges of suspected L1s that are too c
DistToNearestL1    <- nearest(SuspectL1Ranges, L1GRanges)
idxSuspectL1Ranges <- idxSuspectL1Ranges[DistToNearestL1 >= MinDist2L1]
SuspectL1Ranges    <- IslandGRanges[idxSuspectL1Ranges]

# Create a vector of fastq file names
FastQnames <- paste(as.vector(seqnames(SuspectL1Ranges)), 
                    idxSuspectL1Ranges, sep = "_")
FastQPaths <- paste(OutFastQFolder, FastQnames, ".fastq", sep = "")
if (!dir.exists(OutFastQFolder)){
  dir.create(OutFastQFolder)
}

# Write little fastq files per suspected 
cat("*******   Writing little fastq files ...   *******\n")
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
save(list = c("IslandGRanges", "L1GRanges", "ScannedL1Ranges", "ReadsPerL1", 
              "CoverMat", "QuantileMat", "idxRange"), file = OutResults)
