# The following script calculates the average 1000 genome coverage for each bp 
# within reference L1

# Source start script
source('/home/hb54/L1polymORFgit/Scripts/_Start_L1polymORF_AUB.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Run parameters
RunTime   <- '12:00:00'
Mem       <- '50G'
JobsTotal <- 500

# Path to reference data and for 1000 genomes
BedPath_L1  <- '/home/hb54/RefSeqData/L1HSRefRanges_hg38.bed'
Path1000G   <- "/home/hb54/1000GenomeData/"
Files1KG <- ReadVCF("/home/hb54/1000GenomeData/1000genomes.high_coverage.GRCh38DH.alignment.index")

#
Files1KG$X.CRAM <- gsub("ftp:/", "http://", Files1KG$X.CRAM)
RunIDs <- NULL
BamPath1000G <- Files1KG$X.CRAM[2]
for (BamPath1000G in Files1KG$X.CRAM[-1]){
  IndividualID <- strsplit(BamPath1000G, "/")[[1]][10]
  
  cat("*******   Getting coverage for", IndividualID, "    *******\n")
  cat("Job", which(BamPath1000G == Files1KG$X.CRAM),"out of a total of", 
      length(Files1KG$X.CRAM), "\n")
  
  # Define paths
  CoverOutPath   <- paste(Path1000G, "L1Coverage_", IndividualID, sep = "")

  # Construct samtools command to get filtered bam file
  SamToolsCmds <- paste("/home/hb54/samtools-1.11/bin/samtools depth -a -b", BedPath_L1, 
                          BamPath1000G, ">", CoverOutPath)
 

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

save.image(paste(Path1000G, "RunInfo_Coverage_HG38.RData"))
