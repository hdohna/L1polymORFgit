##############################################
#
# General description:
#
#   The following script reads calculates Tajima's D

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
DataFolder     <- "/labs/dflev/hzudohna/1000Genomes/"
ChrLPath           <- '/labs/dflev/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName    <- "/labs/dflev/hzudohna/RefSeqData/L1HS_repeat_table_Hg19.csv"
OutBedPath_L1Neighbor <- paste('/labs/dflev/hzudohna/RefSeqData/L1HS_Neighbor',
                               FlankSize, '_hg19.bed', sep = "")
OutVcfPath         <- '/scratch/users/hzudohna/Variants1000G.vcf'
OutTajimaPrefix    <- paste('/labs/dflev/hzudohna/1000Genomes/L1Neighbor', 
                            FlankSize, sep = "")

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

# Loop through vcf files and estimate Tajima's D 
InFile <- AllFiles[1]
for (InFile in AllFiles){
  InFileSplit <- strsplit(InFile, "\\.")[[1]]
  OutPrefix     <- paste(InFileSplit[1:2], collapse = "_")
  OutPath2    <- paste(c(InFileSplit[2], OutPrefix), collapse = "_")
  cat("Calculate Tajima's D for", InFileSplit[2], "\n\n")
  ScriptName <- paste("L1indelL_Script", InFileSplit[2], sep = "_")
  CreateAndCallSlurmScript(file = paste("runVcftools_L1TajimaD", InFileSplit[2], sep = "_"), 
                           scriptName = paste("Vcftools_L1TajimaD", InFileSplit[2], sep = "_"),
                           SlurmCommandLines = 
                             c("module load vcftools",
                               paste("vcftools --vcf", InFile, 
                                     "--bed", OutBedPath_L1Neighbor,
                                     "--TajimaD", StepSize,
                                     "--out", OutPrefix)),
                           RunTime = '12:00:00',
                           Mem = '200G')
}


# Check whether queue is finished
blnQueueFinished <- CheckQueue(MaxNrTrials = 100,
                       SleepTime   = 30)

# If queue is finished, combine the chromosome-level files into one combined file    
if (blnQueueFinished){
  TajimaFiles <- list.files(DataFolder, pattern =".Tajima.D", 
                         full.names = T)
  TajimaData <- read.delim(TajimaFiles[1])
  for(TajimaFile in TajimaFiles){
    NewData <- read.delim(TajimaFile)
    TajimaData <- rbind(TajimaData, NewData)
  }
  
}

# Create genomic ranges 
DiffBefore <- TajimaData$BIN_START - c(0,TajimaData$BIN_START[-nrow(TajimaData)])
DiffAfter  <- c(TajimaData$BIN_START[-1], 0) - TajimaData$BIN_START
idxStart   <- which(DiffBefore != StepSize)
idxEnd     <- which(DiffAfter != StepSize)
TajimaD_GR <- GRanges(TajimaData, seqnames.field = "CHROM",
                                       start.field = "BIN_START",
                                       end.field = "BIN_START")


TajimaD_GR <- GRanges(seqnames = TajimaData$CHROM[-1], 
                      ranges = IRanges(start = TajimaData$BIN_START))

TajimaOLcount <- countOverlaps(TajimaD_GR, L1Neighborhoods)

save.image('/labs/dflev/hzudohna/1000Genomes/L1_TajimaD.RData')
