# The following script creates objects and commands for calling L1 deletions in
# 1000 genome data using MELT.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Boolean variables for different parts of the workflow
blnRunSam  <- T
blnRunMELT <- T

# Specify run parameters
RunTime <- '12:00:00'
Mem     <- '100G'

# Path to reference data and for 1000 genomes
#RefPath     <- "D:/L1polymORF/Data/"
RefPath      <- "/labs/dflev/hzudohna/RefSeqData/"
RefFilePath  <- "/labs/dflev/hzudohna/RefSeqData/hg19.fa"
Path1000G    <- "/labs/dflev/hzudohna/1000Genomes/"
Path1000GSim <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Create bed file of ranges around each reference L1
ChrNames <- as.vector(seqnames(L1GRanges))
ChrNrs   <- substr(ChrNames, 4, nchar(ChrNames))
L1NeighborRanges <- GRanges(seqnames = ChrNrs, 
                            IRanges(start = start(L1GRanges) - 500,
                                    end   = end(L1GRanges) + 500))
L1NeighborbedPath <- paste(RefPath, "L1NeighborRanges.bed", sep = "")
L1bedPath         <- paste(RefPath, "L1Ranges.bed", sep = "")
export.bed(L1NeighborRanges, L1NeighborbedPath) 

L1Ranges_shortName <- GRanges(seqnames = ChrNrs, 
                              IRanges(start = start(L1GRanges),
                                      end   = end(L1GRanges)))
export.bed(L1Ranges_shortName, L1bedPath) 


# Create bed file of ranges around each 1000 genome L1
ChrNames <- as.vector(seqnames(L1_1000G_GR_hg19))
L1NeighborRanges_1000G <- GRanges(seqnames = ChrNames, 
                            IRanges(start = start(L1_1000G_GR_hg19) - 500,
                                    end   = end(L1_1000G_GR_hg19) + 500))
L1NeighborbedPath_1000G <- paste(RefPath, "L1NeighborRanges_1000G_chr.bed", sep = "")
L1bedPath_1000G         <- paste(RefPath, "L1Ranges_1000G.bed", sep = "")
#export.bed(L1NeighborRanges_1000G, L1NeighborbedPath_1000G) 

export.bed(L1NeighborRanges_1000G, L1NeighborbedPath_1000G) 

# Get paths to bam files
BamFiles <- list.files(Path1000GSim, pattern = "_sorted.bam", full.names = T)
BamFiles <- BamFiles[-grep("_sorted.bam.", BamFiles)]
BamFiles_Filtered <- list.files(Path1000GSim, pattern = "Filtered.bam", full.names = T)
BamFiles_Filtered <- BamFiles_Filtered[-grep("Filtered.bam.", BamFiles_Filtered)]
BamFiles_Filtered

########################################################################################
# Run MELT for low coverage bam files
BamFile <- BamFiles[1]
for (BamFile in BamFiles){
  
  # Get ID
  Split1 <- strsplit(BamFile, "_")[[1]]

  # Create folder for MELT
  NewDir <- paste(Path1000GSim, Split1[4], "_Filtered", sep = "")
  
  # Define paths
  BamOutPath   <- gsub("_sorted.bam", "_sortedL1Filtered.bam", BamFile)
  
  # Construct samtools command to get filtered bam file
  SamToolsCmds <- c("module load samtools",
                    paste("samtools view -b -h", BamFile, "-L", 
                          L1NeighborbedPath_1000G, "| samtools sort -o",
                          BamOutPath),
                    paste("samtools index", BamOutPath))
  
  # Construct MELT command
  MkDirCmd <- NULL
  if (!dir.exists(NewDir)){
    paste("mkdir", NewDir)
  }
  MELTCmds <- c(MkDirCmd,
                "module load bowtie2",
                paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
                      BamOutPath, 
                      "-c 7",
                      "-h", RefFilePath,
                      "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                      "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                      "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                      "-w", NewDir))
  
  
  # Create a list of all commands
  # Create a list of all commands
  CmdList <- list(SamToolsCmds, MELTCmds)
  AllCmds <- unlist(CmdList[c(blnRunSam, blnRunMELT)])

  # Create script name and run script
  ScriptName <- paste("Sim_Bwa_MELTScript", Split1[4], "Filter", sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           RunTime = RunTime,
                           Mem = Mem,
                           SlurmCommandLines = AllCmds)
}

