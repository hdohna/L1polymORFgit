# The following script obtains barcodes of reads aligning to flanks of catalgue
# L1 (produced by script 'extract10Xbarcodes.R') and counts the number of reads 
# that align to L1 


# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(ShortRead)
library(csaw)
library(rtracklayer)

# Specify bam file path
InBamfilePath <- "/share/diskarray3/hzudohna/10XData/NA12878_WGS_phased_possorted.bam"
OutFilePrefix <- "/share/diskarray3/hzudohna/10XData/L1_"
BamFolder     <- "/share/diskarray3/hzudohna/10XData/"
BamPrefix     <- "L1_"
L1HSConsensus <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
L1HSCatRef   <- "/home/hzudohna/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fas"

# Specify commands
AlignCommand <- '/home/txw/bwa/bwa-0.7.12/bwa mem -t6'

# Read in table with known L1 
L1Catalogue <- read.csv("/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv", 
                        as.is = T)

# Get bam files with specified prefix, loop over files and extract accession numbers
BamFiles <- list.files(BamFolder, pattern = BamPrefix, full.names = T)
BamFiles <- BamFiles[grep(".bam", BamFiles)]
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}
if (length(grep("_aln2L1", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep("_aln2L1", BamFiles)]
}
BamFile <- BamFiles[2]

# Initialize lists that keep track of reads
ReadListL1Consens   <- list()
ReadListL1Catalogue <- list()

# Create genomic ranges for consensus and catalogue L1
GRL1consens <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
              ranges = IRanges(start = 1, end = 6000))
L1CatLines <- readLines(L1HSCatRef)
L1CatNames <- grep(">", L1CatLines, value = T)
L1CatNames <- gsub(">", "", L1CatNames)
GRL1Catalogue <- GRanges(seqnames = L1CatNames, 
                       ranges = IRanges(start = 1, end = 16000))

# Map reads from bam file to catalogue and consensus L1
for (BamFile in BamFiles){
  
  cat("*********  Processing", BamFile, "\n\n")
  
  # Create file names
  FastqFile <- gsub(".bam", ".fastq", BamFile)
  SamFileL1Consens <- paste(strsplit(BamFile, ".bam")[[1]][1], 
                            "_aln2L1Consens.sam", sep = "")
  SamFileL1Catalogue <- paste(strsplit(BamFile, ".bam")[[1]][1], 
                            "_aln2L1Catalogue.sam", sep = "")
  BamFileL1Consens   <- gsub(".sam", ".bam", SamFileL1Consens)
  BamFileL1Catalogue <- gsub(".sam", ".bam", SamFileL1Catalogue)
  BamDestL1Consens   <- gsub(".sam", "", SamFileL1Consens)
  BamDestL1Catalogue <- gsub(".sam", "", SamFileL1Catalogue)
  BamFileL1CatFilter <- gsub(".bam", "Filtered.bam", BamFileL1Catalogue)
  
  # Turn bam file into fastq file
  Cmd <- paste("/home/txw/samtools/samtools-1.2/samtools bam2fq", BamFile,
               ">", FastqFile)
  system(Cmd)
  
  # Align fastq file to L1 consensus
  AlnCmd <- paste(AlignCommand, L1HSConsensus, FastqFile, ">", SamFileL1Consens)
  system(AlnCmd)
  
  # Align fastq file to L1 catalogue
  AlnCmd <- paste(AlignCommand, L1HSCatRef, FastqFile, ">", SamFileL1Catalogue)
  system(AlnCmd)
  
  # Turn sam files into bam files
  asBam(SamFileL1Consens,   BamDestL1Consens)
  asBam(SamFileL1Catalogue, BamDestL1Catalogue, overwrite = T)
  
  # Extract reads aligned to L1HS consensus 
  ReadListL1Consens   <- c(ReadListL1Consens, 
                           extractReads(BamFileL1Consens, GRL1consens))
  
  # Filter and extract reads aligned to L1HS L1 catalogue
  paramFilter  <- ScanBamParam(mapqFilter = 1)
  filterBam(BamFileL1Catalogue, BamFileL1CatFilter, param = paramFilter)
  RLC <- lapply(GRL1Catalogue, function(x)extractReads(BamFileL1CatFilter, x))
  ReadListL1Catalogue <- c(ReadListL1Catalogue, RLC)
}

# Save lists with reads
cat("*********   Saving read lists    **********")
save(list = c("ReadListL1Consens", "ReadListL1Catalogue"),
     file = "/share/diskarray3/hzudohna/10XData/L1Catalogue10XReadLists.RData")
