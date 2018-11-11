# The following script creates objects and commands for calling L1 deletions in
# 1000 genome data using MELT.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Path to reference data and for 1000 genomes
#RefPath   <- "D:/L1polymORF/Data/"
RefPath   <- "/labs/dflev/hzudohna/RefSeqData/"
Path1000G <- "/labs/dflev/hzudohna/1000Genomes/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Create bed file of ranges around each L1
ChrNames <- as.vector(seqnames(L1GRanges))
ChrNrs   <- substr(ChrNames, 4, nchar(ChrNames))
L1NeighborRanges <- GRanges(seqnames = ChrNrs, 
                            IRanges(start = start(L1GRanges) - 5000,
                                    end   = end(L1GRanges) + 5000))
L1NeighborbedPath <- paste(RefPath, "L1NeighborRanges.bed", sep = "")
L1bedPath         <- paste(RefPath, "L1Ranges.bed", sep = "")
export.bed(L1NeighborRanges, L1NeighborbedPath) 

L1Ranges_shortName <- GRanges(seqnames = ChrNrs, 
                              IRanges(start = start(L1GRanges),
                                      end   = end(L1GRanges)))
export.bed(L1Ranges_shortName, L1bedPath) 

# Specify the general path to 1000 genome bam file
BamPath1000G_General <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/IndividualID/alignment/IndividualID.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
DirPath1000G_General <- "ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/IndividualID/alignment/"

# Get all files in the 1000 genome path
Files1000G <- list.files(Path1000G)
blnNotAnalyzed <- sapply(SampleColumns, function(x) length(grep(x, Files1000G)) == 0)
IndividualID <- SampleColumns[blnNotAnalyzed][3]
SampleColumns[blnNotAnalyzed]

for (IndividualID in SampleColumns[blnNotAnalyzed]){
  
  # Define paths
  DirPath1000G <- gsub("IndividualID", IndividualID, DirPath1000G_General)
  BamOutPath   <- paste(Path1000G, "L1Filtered_", IndividualID, ".bam", sep = "")
  
  # Create curl command to get file names on ftp directories
  DirListPath  <- paste(Path1000G, "FileList.txt", sep = "")
  CurlCmd <- paste("curl", paste("ftp://", DirPath1000G, sep = ""),
                   ">", DirListPath)
  system(CurlCmd)

  # Get name of bam file 
  DirListLines <- readLines(DirListPath)
  FileNames <- sapply(DirListLines, function(x){
    SplitLines <- strsplit(x, " ")[[1]]
    SplitLines[length(SplitLines)]
  })
  FileNames   <- FileNames[grep("\\.mapped", FileNames)]
  BamFileName <- FileNames[-grep("bam.b", FileNames)]
  BamPath1000G <- paste("ftp://", DirPath1000G, BamFileName, sep = "")

  # Construct samtools command to get filtered bam file
  SamToolsCmds <- c("module load samtools",
                    paste("samtools view -b -h", BamPath1000G, "-L", 
                          L1NeighborbedPath, "| samtools sort -o",
                          BamOutPath),
                    paste("samtools index", BamOutPath))
  
  # Construct MELT command
  MELTCmds <- paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Deletion-Genotype -bamfile",
                    BamOutPath,"-h /labs/dflev/hzudohna/RefSeqData/hg19masked.fa",
                    "-bed", L1bedPath,
                    "-w /labs/dflev/hzudohna/1000Genomes/")
  
  # Command to remove bam files
  RemCmds <- c(paste("rm", BamOutPath), paste("rm ", BamOutPath, ".bai", sep = ""))
  
  ScriptName <- paste("ME_DEL_Script", IndividualID, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           SlurmCommandLines = c(SamToolsCmds, 
                                                 "echo 'bam file retrieved'", 
                                                 MELTCmds, 
                                                 "echo 'MELT completed'",
                                                 RemCmds,
                                                 "echo 'bam files removed'"))
  
}



