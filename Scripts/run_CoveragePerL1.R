# The following script calculates the average 1000 genome coverage for each bp 
# within reference L1

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Run parameters
RunTime   <- '12:00:00'
Mem       <- '100G'
JobsTotal <- 500

# Path to reference data and for 1000 genomes
RefPath     <- "/labs/dflev/hzudohna/RefSeqData/"
DataFolder         <- "/labs/dflev/hzudohna/1000Genomes/"
ChrLPath           <- '/labs/dflev/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName    <- "/labs/dflev/hzudohna/RefSeqData/L1HS_repeat_table_Hg19.csv"
OutBedPath_L1      <- '/labs/dflev/hzudohna/RefSeqData/L1HS_hg19.bed'

RefFilePath <- "/reference/RefGenomes/1000genomes/hs37d5/hs37d5.fa"
Path1000G   <- "/labs/dflev/hzudohna/1000Genomes/"

# Path for scratch storage
ScratchPath <- "/scratch/users/hzudohna/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName, as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))

# Create GRanges objects with L1 Seqences
L1GRanges <- makeGRangesFromDataFrame(L1Table, seqnames.field = "ChrNr",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")

# Specify the general path to 1000 genome bam file
BamPath1000G_General <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/IndividualID.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
DirPath1000G_General <- "ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/IndividualID/alignment/"

#######################################
#                                     #
#    Export data and run samtools     #
#                                     #
#######################################

# Export bed ranges of L1
export.bed(L1GRanges, con = OutBedPath_L1)

# Get names of potential coverage files
CoverFiles <- paste(Path1000G, "L1Coverage_", SampleColumns, sep = "")

RunIDs <- NULL
IndividualID <- SampleColumns[1]
IDs2Analyze <- SampleColumns[!file.exists(CoverFiles)]
for (IndividualID in IDs2Analyze){
  cat("*******   Getting coverage for", IndividualID, "    *******\n")
  cat("Job", which(IDs2Analyze == IndividualID),"out of a total of", 
      length(IDs2Analyze), "\n")
  
  # Define paths
  DirPath1000G <- gsub("IndividualID", IndividualID, DirPath1000G_General)
  CoverOutPath   <- paste(Path1000G, "L1Coverage_", IndividualID, sep = "")

  # Create curl command to get file names on ftp directories
  DirListPath  <- paste(Path1000G, IndividualID, "_FileList.txt", sep = "")
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
  BamFileName <- FileNames[-grep("bam.", FileNames)]
  BaiFileName <- FileNames[grep(".bai", FileNames)]
  BamPath1000G <- paste("ftp://", DirPath1000G, BamFileName, sep = "")
  BaiPath1000G <- paste("ftp://", DirPath1000G, BaiFileName, sep = "")
  
  # Construct samtools command to get filtered bam file
  SamToolsCmds <- c("module load samtools",
                    paste("samtools depth -a -b", OutBedPath_L1, 
                          BamPath1000G, ">", CoverOutPath))
 

  # Create script name and launch script
  ScriptName <- paste("L1Coverage_", IndividualID, sep = "_")
  RunList <- CreateAndCallSlurmScript(file = ScriptName, 
                           scriptName = IndividualID,
                           SlurmCommandLines = SamToolsCmds,
                           RunTime = RunTime,
                           Mem = Mem)
  cat(RunList$RunMessage, "\n")
  RunIDs <- c(RunIDs, RunList$RunID)
}
