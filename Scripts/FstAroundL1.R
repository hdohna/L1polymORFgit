##############################################
#
# General description:
#
#   The following script reads calculates Fst around L1s

# Input:
#

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
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Specify flanking region and step size
StepSize   <- 10000
FlankSize  <- 500000

# Specify paths
DataFolder         <- "D:/L1polymORF/Data/"
DataFolder         <- "/labs/dflev/hzudohna/1000Genomes/"
ChrLPath           <- '/labs/dflev/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName    <- "/labs/dflev/hzudohna/RefSeqData/L1HS_repeat_table_Hg19.csv"
OutBedPath_L1Neighbor <- paste('/labs/dflev/hzudohna/RefSeqData/L1HS_Neighbor',
                               FlankSize, '_hg19.bed', sep = "")
OutVcfPath         <- '/scratch/users/hzudohna/Variants1000G.vcf'
OutFstPrefix    <- paste('/labs/dflev/hzudohna/1000Genomes/L1Neighbor', 
                            FlankSize, sep = "")
PopPrefix <- "1000G_Pop_"

# Get names of all sample population files 
SampleInfo        <- read.delim("D:/L1polymORF/Data/1000GenomeSampleInfo.txt")
SamplePopFiles    <- list.files(DataFolder, pattern = PopPrefix, full.names = T)
SamplePopFileInfo <- file.info(SamplePopFiles)
SamplePopFiles    <- SamplePopFiles[SamplePopFileInfo$size > 0]

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Load vector with chromosome lengths
load(ChrLPath)

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName, as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))

# Create GRanges objects with L1 Seqences
L1GRanges <- makeGRangesFromDataFrame(L1Table, seqnames.field = "ChrNr",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")

# Resize L1 GRanges to get neighborhoods
L1Neighborhoods <- resize(L1GRanges, FlankSize, fix = "center")

#######################################
#                                     #
#    Export data and run vcftools     #
#                                     #
#######################################

# Export bed ranges
export.bed(L1Neighborhoods, con = OutBedPath_L1Neighbor)

# Create a vector of fst commands per population
FstCmds <- paste("--weir-fst-pop", SamplePopFiles, collapse = " ")

# Loop through chromosovcf files and estimate Fst for neighborhoods around L1
InFile <- AllFiles[1]
RunOverview <- data.frame()
cat("\n**********   Fst for neighborhoods arond L1   *****\n")
for (InFile in AllFiles){
  
  # Create and submit a script that calculates Fst for all populations
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutPrefix   <- paste(c(InFileSplit[1:2], "AllPops"), collapse = "_")
  OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
  cat("Calculate Fst for", InFileSplit[2], "\n\n")
  ScriptName <- paste("L1indelL_Script", InFileSplit[2], sep = "_")
  RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1Fst", InFileSplit[2], sep = "_"), 
                           scriptName = paste("Vcftools_L1Fst", InFileSplit[2], sep = "_"),
                           SlurmCommandLines = 
                             c("module load vcftools",
                               paste("vcftools --vcf", InFile, 
                                     "--bed", OutBedPath_L1Neighbor,
                                     FstCmds,
                                     "--fst-window-step", StepSize,
                                     "--out", OutPrefix)),
                           RunTime = '12:00:00',
                           Mem = '200G')$RunID
  NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                        Pop = "AllPops", Range = "L1Neighborhood")
  RunOverview <- rbind(RunOverview, NewData)
  Sys.sleep(3)
}

# Loop through vcf files and estimate Fst 
InFile <- AllFiles[1]
RunOverview <- data.frame()
cat("\n**********   Fst for L1 regions    *****\n")
for (InFile in AllFiles){
  
  # Create and submit a script that calculates Fst for all populations
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutPrefix   <- paste(c(InFileSplit[1:2], "AllPops"), collapse = "_")
  OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
  cat("Calculate Fst for", InFileSplit[2], "\n\n")
  ScriptName <- paste("L1indelL_Script", InFileSplit[2], sep = "_")
  RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1Fst", InFileSplit[2], sep = "_"), 
                                    scriptName = paste("Vcftools_L1Fst", InFileSplit[2], sep = "_"),
                                    SlurmCommandLines = 
                                      c("module load vcftools",
                                        paste("vcftools --vcf", InFile, 
                                              "--bed", OutBedPath_L1Neighbor,
                                              FstCmds,
                                              "--fst-window-step", StepSize,
                                              "--out", OutPrefix)),
                                    RunTime = '12:00:00',
                                    Mem = '200G')$RunID
  NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                        Pop = "AllPops", Range = "L1Neighborhood")
  RunOverview <- rbind(RunOverview, NewData)
  Sys.sleep(3)
}

# Check whether queue is finished
blnQueueFinished <- CheckQueue(MaxNrTrials = 200,
                       SleepTime   = 60,
                       JobIDs = RunOverview$RunID[!is.na(RunOverview$RunID)])

# Add column about job completion once jobs are finished
if (blnQueueFinished){
  RunOverview$blnCompleted <- CheckJobCompletion(JobIDs = RunOverview$RunID)
}


# If queue is finished, combine the chromosome-level files into one combined file    
if (blnQueueFinished){
  FstFiles <- list.files(DataFolder, pattern =".weir.fst", 
                         full.names = T)
  FstData <- read.delim(FstFiles[1])
  for(FstFile in FstFiles){
    NewData     <- read.delim(FstFile)
    Split1      <- strsplit(FstFile, "\\/")[[1]]
    FilePart    <- Split1[length(Split1)]
    Split2      <- strsplit(FilePart, "\\.")[[1]]
    Split3      <- strsplit(Split2[1], "\\_")[[1]]
    NewData$Pop <- Split3[length(Split3)]
    FstData  <- rbind(FstData, NewData)
  }
  
}

# Create genomic ranges 
DiffBefore <- FstData$BIN_START - c(0,FstData$BIN_START[-nrow(FstData)])
DiffAfter  <- c(FstData$BIN_START[-1], 0) - FstData$BIN_START
idxStart   <- which(DiffBefore != StepSize)
idxEnd     <- which(DiffAfter != StepSize)
Fst_GR <- GRanges(FstData, seqnames.field = "CHROM",
                                       start.field = "BIN_START",
                                       end.field = "BIN_START")


Fst_GR <- GRanges(seqnames = FstData$CHROM[-1], 
                      ranges = IRanges(start = FstData$BIN_START))

FstOLcount <- countOverlaps(Fst_GR, L1Neighborhoods)

save.image('/labs/dflev/hzudohna/1000Genomes/L1_Fst.RData')
