# The following script creates a fasta file that contain the first 2000 bp of
# L1s in reference and non-reference insertions

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)

# Load ranges
load("D:/L1polymORF/Data/BZ_L1Ranges.RData")

# Read L1 consensus sequence

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv", as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1),]

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
                                  ChainFilePath = "D:/L1polymORF/Data/hg38ToHg19.over.chain")
L1CatalogGR <- LiftOverList$GRCatalogue_hg19# Specify folder 
OutFolderName_NonRef <- "D:/L1polymORF/Data/BZ_NonRef"

# get names of newly created bam files
FileNames <- list.files(OutFolderName_NonRef, pattern = ".bam",
                        full.names = T)
FileNames <- FileNames[-grep(".bam.", FileNames)]

# Loop through file names and read in bam files of reads mapped to L1
ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))

# Count the number of reads mapped
NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
  sum(!is.na(x[[1]]$pos))
})
idxFilesWithReads <- which(NrMapped2L1 > 0)
FilesWithReads <- FileNames[idxFilesWithReads]

# Get read list per peak
ReadListPerPeak <- lapply(idxFilesWithReads, function(x) {
  ScannedL1Ranges[[x]][[1]]
})

RL <- ReadListPerPeak[[2]]
i <- 1
RL$pos
SeqList <- lapply(1:3, function(i) {
  StartSeq <- rep("-", RL$pos[i])
  RSeq <- SeqFromCigar(RL$cigar[i], RL$seq[i])
  c(StartSeq, RSeq)
  })
SeqLMaxLength <- max(sapply(SeqList, length))
SeqMat <- t(sapply(SeqList, function (x) c(x, rep("-", SeqLMaxLength - length(x)))))

