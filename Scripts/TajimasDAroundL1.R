##############################################
#
# General description:
#
#   The following script reads calculates Tajima's D around L1s

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
OutBedPath_L1      <- '/labs/dflev/hzudohna/RefSeqData/L1HS_hg19.bed'
OutBedPath_L1Neighbor <- paste('/labs/dflev/hzudohna/RefSeqData/L1HS_Neighbor',
                               FlankSize, '_hg19.bed', sep = "")
OutVcfPath         <- '/scratch/users/hzudohna/Variants1000G.vcf'
OutTajimaPrefix    <- paste('/labs/dflev/hzudohna/1000Genomes/L1Neighbor', 
                            FlankSize, sep = "")
PopPrefix <- "1000G_Pop_"

# Get names of all sample population files 
SampleInfo        <- read.delim("/labs/dflev/hzudohna/1000Genomes/1000GenomeSampleInfo.txt")
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
export.bed(L1GRanges, con = OutBedPath_L1)
export.bed(L1Neighborhoods, con = OutBedPath_L1Neighbor)

# Loop through chromosomes and estimate Tajima's D for neighborhoods and all 
# populations
RunOverview <- data.frame()
cat("\n**********   Tajima's D for neighborhoods and all populations    *****\n")
for (InFile in AllFiles){
  
  # Create and submit a script that calculates Tajima's D for all populations
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutPrefix   <- paste(c(InFileSplit[1:2], FlankSize, "AllPops"), collapse = "_")
  OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
  cat("Calculate Tajima's D for", InFileSplit[2], "\n")
  ScriptName <- paste("Vcftools_L1TajimaD", InFileSplit[2], FlankSize, "AllPops", sep = "_")
  RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1TajimaD", InFileSplit[2], sep = "_"), 
                           scriptName = ScriptName,
                           SlurmCommandLines = 
                             c("module load vcftools",
                               paste("vcftools --vcf", InFile, 
                                     "--bed", OutBedPath_L1Neighbor,
                                     "--TajimaD", StepSize,
                                     "--out", OutPrefix)),
                           RunTime = '12:00:00',
                           Mem = '200G')$RunID
  NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                        Pop = "AllPops", Range = "L1Neighborhood")
  RunOverview <- rbind(RunOverview, NewData)
  Sys.sleep(3)
}


# Loop through chromosomes and estimate Tajima's D for L1 only and all 
# populations
cat("\n**********   Tajima's D for L1 and all populations    *****\n")
for (InFile in AllFiles){
  
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  cat("Calculate Tajima's D for", InFileSplit[2], "\n")
  OutPrefix   <- paste(c(InFileSplit[1:2],  "L1_AllPops"), collapse = "_")
  ScriptName <- paste("Vcftools_L1TajimaD", InFileSplit[2], "L1_AllPops", sep = "_")
  RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1TajimaD", InFileSplit[2], "L1", sep = "_"), 
                           scriptName = ScriptName,
                           SlurmCommandLines = 
                             c("module load vcftools",
                               paste("vcftools --vcf", InFile, 
                                     "--bed", OutBedPath_L1,
                                     "--TajimaD", 10000,
                                     "--out", OutPrefix)),
                           RunTime = '12:00:00',
                           Mem = '200G')$RunID
  NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                        Pop = "AllPops", Range = "L1")
  RunOverview <- rbind(RunOverview, NewData)
  Sys.sleep(3)
  
}  

# Loop through chromosomes and estimate Tajima's D for neighborhoods and each 
# population 
cat("\n**********   Tajima's D for neighborhoods and individual populations    *****\n")
for (InFile in AllFiles){
  
  # Get parts of input file
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  cat("Calculate Tajima's D for", InFileSplit[2], "\n")
  
  # Loop over populations and estimate Tajima's D for each population
  for (PopFile in SamplePopFiles){
    Pop         <- strsplit(PopFile, PopPrefix)[[1]][2]
    OutPrefix   <- paste(c(InFileSplit[1:2], FlankSize, Pop), collapse = "_")
    OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
    
    cat("Calculate Tajima's D for", InFileSplit[2], Pop, "\n")
    ScriptName <- paste("Vcftools_L1TajimaD", InFileSplit[2], FlankSize, Pop, sep = "_")
    RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1TajimaD", InFileSplit[2], FlankSize, Pop, sep = "_"), 
                             scriptName = ScriptName,
                             SlurmCommandLines = 
                               c("module load vcftools",
                                 paste("vcftools --vcf", InFile, 
                                       "--bed", OutBedPath_L1Neighbor,
                                       "--keep", PopFile,
                                       "--TajimaD", StepSize,
                                       "--out", OutPrefix)),
                             RunTime = '12:00:00',
                             Mem = '200G')$RunID
    NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                          Pop = Pop, Range = "L1Neighborhood")
    RunOverview <- rbind(RunOverview, NewData)
    Sys.sleep(3)
    
  }
}

# Loop through chromosomes and estimate Tajima's D for L1s and each 
# population 
cat("\n**********   Tajima's D for L1 and individual populations    *****\n")
for (InFile in AllFiles){
  
  # Get parts of input file
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  cat("Calculate Tajima's D for", InFileSplit[2], "\n")
  
  # Loop over populations and estimate Tajima's D for each population
  for (PopFile in SamplePopFiles){
    Pop         <- strsplit(PopFile, PopPrefix)[[1]][2]
    cat("Calculate Tajima's D for", InFileSplit[2], Pop, "\n")
    OutPrefix   <- paste(c(InFileSplit[1:2], "L1", Pop), collapse = "_")
    ScriptName <- paste("Vcftools_L1TajimaD", InFileSplit[2], "L1", Pop, sep = "_")
    RunID <- CreateAndCallSlurmScript(file = paste("runVcftools_L1TajimaD", InFileSplit[2],"L1", Pop, 
                                          sep = "_"), 
                             scriptName = ScriptName,
                             SlurmCommandLines = 
                               c("module load vcftools",
                                 paste("vcftools --vcf", InFile, 
                                       "--bed", OutBedPath_L1,
                                       "--keep", PopFile,
                                       "--TajimaD", 10000,
                                       "--out", OutPrefix)),
                             RunTime = '12:00:00',
                             Mem = '200G')$RunID
    NewData <- data.frame(RunID = RunID, Chrom = InFileSplit[2],
                          Pop = Pop, Range = "L1")
    RunOverview <- rbind(RunOverview, NewData)
    Sys.sleep(1)
    
  }
}


# Check whether queue is finished
blnQueueFinished <- CheckQueue(MaxNrTrials = 300,
                       SleepTime   = 60,
                       JobIDs = RunOverview$RunID[!is.na(RunOverview$RunID)])

# Add column about job completion once jobs are finished
# if (blnQueueFinished){
#   RunOverview$blnCompleted <- CheckJobCompletion(JobIDs = RunOverview$RunID)
# }
  
# If queue is finished, combine the chromosome-level files into one combined file    
if (blnQueueFinished){
  TajimaFiles <- list.files(DataFolder, pattern =".Tajima.D", 
                         full.names = T)
  TajimaData <- read.delim(TajimaFiles[1])
  Split1      <- strsplit(TajimaFile[1], "\\/")[[1]]
  FilePart    <- Split1[length(Split1)]
  Split2      <- strsplit(FilePart, "\\.")[[1]]
  Split3      <- strsplit(Split2[1], "\\_")[[1]]
  TajimaData$Pop <- Split3[length(Split3)]
  for(TajimaFile in TajimaFiles){
    NewData     <- read.delim(TajimaFile)
    Split1      <- strsplit(TajimaFile, "\\/")[[1]]
    FilePart    <- Split1[length(Split1)]
    Split2      <- strsplit(FilePart, "\\.")[[1]]
    Split3      <- strsplit(Split2[1], "\\_")[[1]]
    NewData$Pop <- Split3[length(Split3)]
    TajimaData  <- rbind(TajimaData, NewData)
  }
  
}

save(list = "TajimaData", file = '/labs/dflev/hzudohna/1000Genomes/L1_TajimaData.RData')

# Create genomic ranges 
# DiffBefore <- TajimaData$BIN_START - c(0,TajimaData$BIN_START[-nrow(TajimaData)])
# DiffAfter  <- c(TajimaData$BIN_START[-1], 0) - TajimaData$BIN_START
# idxStart   <- which(DiffBefore != StepSize)
# idxEnd     <- which(DiffAfter != StepSize)
# TajimaD_GR <- GRanges(TajimaData, seqnames.field = "CHROM",
#                                        start.field = "BIN_START",
#                                        end.field = "BIN_START")
# 
# 
# TajimaD_GR <- GRanges(seqnames = TajimaData$CHROM[-1], 
#                       ranges = IRanges(start = TajimaData$BIN_START))
# 
# TajimaOLcount <- countOverlaps(TajimaD_GR, L1Neighborhoods)

save.image('/labs/dflev/hzudohna/1000Genomes/L1_TajimaD.RData')
