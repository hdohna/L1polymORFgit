##############################################
#
# General description:
#
#   The following script reads a bam file, ranges of known L1 on the 
#   reference genome 

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
CoverComparePlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverComparison.pdf'
OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference.Rdata'

# Suffices for alignment files created by BWA
SamSuffix <- "_aln.sam"
BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")

# BWA command (options can be added here)
BWAcommand <- 'bwa mem'

# Peak calling parameters
MinMaxCover <- 40    # minimum maximum coverage to be called a peak 
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
#    Define ranges of known and suspected      #
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


# Save results
cat("*******  Saving results ...   *******\n")
save(list = c("IslandGRanges", "L1GRanges", "ScannedL1Ranges", "ReadsPerL1", 
              "CoverMat", "QuantileMat", "idxRange"), file = OutResults)
