# The script below preprocesses 1000 genome data for analysis with selscan
# Data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

########################################
#                                      #
#         Set parameters               #
#                                      #
########################################

# Specify parameters
WindowWidth  <- 10^5
MinNrCarrier <- 5
NrInfoCol    <- 9
MaxNrTrials  <- 20
SleepTime    <- 30
ThinDist     <- 1
NrJobsPerBatch <- 50
WaitBetwJobs <- 100
DataFolder   <- "/labs/dflev/hzudohna/1000Genomes/"
FilePattern  <- "genotypes.vcf"
BedPath      <- "/labs/dflev/hzudohna/1000Genomes/L1WindowSubset.bed"


########################################
#                                      #
#       Process L1 data                #
#                                      #
########################################

# Read in 1000 genome L1 table
L1_1000G <- read.delim("/srv/gsfs0/projects/levinson/hzudohna/1000Genomes/L1_1000G_withGenoNum",
                       header = T, sep = " ")

# Subset L1 file to retain only entries with enough carriers
blnL1 <- L1_1000G[,(NrInfoCol + 1):ncol(L1_1000G)] > 0
blnEnough <- rowSums(blnL1, na.rm = T) >= MinNrCarrier
sum(blnEnough)

# Create a bed file with windows around each L1
L1WindowsBed <- data.frame(chrom  = L1_1000G$CHROM[blnEnough],
                           chromStart = L1_1000G$POS[blnEnough] - WindowWidth,
                           chromEnd   = L1_1000G$POS[blnEnough] + WindowWidth
                           )

# Get index of overlapping rows
NR    <- nrow(L1WindowsBed)
idxOL <- which((L1WindowsBed$chromStart[-1] <= L1WindowsBed$chromEnd[-NR]) &
  (L1WindowsBed$chrom[-1] == L1WindowsBed$chrom[-NR]))

# Remove overlapping rows
L1WindowsBed$chromEnd[idxOL] <- L1WindowsBed$chromEnd[idxOL + 1]
L1WindowsBed <- L1WindowsBed[-(idxOL + 1), ]

# Write out bed file
write.table(L1WindowsBed, BedPath, quote = F, row.names = F)

# Create a map table
MapTable <- cbind(L1_1000G[blnEnough, c("CHROM", "ID")], GenPos = 0.01, 
                  L1_1000G[blnEnough, c("POS")])

# Write out an example map file
MapFile <- paste(DataFolder, "SelscanMap", sep = "")
write.table(MapTable, file = MapFile,
            row.names = F, col.names = F, quote = F, sep = " ")

########################################
#                                      #
#      Subset vcf around L1            #
#                                      #
########################################

# Get file names, loop over files and do the filtering
AllFiles <- list.files(DataFolder, pattern = FilePattern, full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
AllFiles

# Loop over file names, read file and append to existing
cat("\n*********    Subsetting vcf file around each L1  ************\n")
VcfFile <- AllFiles[1]

for (VcfFile in AllFiles){
  cat("Processing", VcfFile, "\n")
  Chrom <- strsplit(VcfFile, "\\.")[[1]][2]
  OutFile <- paste(DataFolder, Chrom, "_L1Windowsubset", sep = "")
  VcfCmd <- c("module load vcftools", 
              paste("vcftools --vcf", VcfFile, "--bed", BedPath, "--recode",
                    "--thin", ThinDist,
              "--out", OutFile))
  ScriptName <- paste("L1Window", Chrom, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName,
                           SlurmHeaderLines = c('#!/bin/bash', 
                                                '#SBATCH --account=dflev', 
                                                '#SBATCH --time=1:00:00', 
                                                '#SBATCH --job-name="Selscan"', 
                                                '#SBATCH --nodes=1', 
                                                '#SBATCH --ntasks=1', 
                                                '#SBATCH --cpus-per-task=1', 
                                                '#SBATCH --mem=100G' 
                           ),
                           SlurmCommandLines = VcfCmd, 
                           scriptName = ScriptName) 
}

# Check whether queue is done (requires that no other batch jobs are running)
QueueFinished <- CheckQueue(MaxNrTrials = MaxNrTrials,
                            SleepTime   = SleepTime)
if (!QueueFinished){
  stop("Subsetting vcf around L1 did not finish within alotted time!")
}

########################################
#                                      #
#     Remove multi-allelic variants    #
#                                      #
########################################

cat("\n*********  Removing multi-allelic variants  ************\n")

# Get file names, loop over files and do the filtering
SubsetFiles <- list.files(DataFolder, pattern = "_L1Windowsubset.recode.vcf", 
                          full.names = T)

# Loop over file names, read file and append to existing
InFile  <- SubsetFiles[1]

for (InFile in SubsetFiles){
  cat("Processing", InFile, "\n")
  Chrom   <- strsplit(strsplit(InFile, "\\_")[[1]][1], "\\//")[[1]][2]
  OutFile <- gsub(".recode", "NoMulti", InFile)
  GrepCmd <- paste('grep -e MULTI_ALLELIC -e \'2|\' -e \'|2\' -v', InFile, ">",
                   OutFile)
  ScriptName <- paste("grepScript", Chrom, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName,
                           SlurmHeaderLines = c('#!/bin/bash', 
                                                '#SBATCH --account=dflev', 
                                                '#SBATCH --time=1:00:00', 
                                                '#SBATCH --job-name="Selscan"', 
                                                '#SBATCH --nodes=1', 
                                                '#SBATCH --ntasks=1', 
                                                '#SBATCH --cpus-per-task=1', 
                                                '#SBATCH --mem=100G' 
                           ),
                           SlurmCommandLines = GrepCmd, 
                           scriptName = ScriptName) 
}

# Check whether queue is done 
QueueFinished <- CheckQueue(MaxNrTrials = MaxNrTrials,
                            SleepTime   = SleepTime)
if (!QueueFinished){
  stop("Removing multi-allelic variants did not finish within alotted time!")
}

