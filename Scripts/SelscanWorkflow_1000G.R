# The script below runs selscan on 1000 Genome data
# It uses pre-processed files that were created by the script `SelscanPreprocessing_1000G`
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
MaxNrTrials  <- 30
SleepTime    <- 10
ThinDist     <- 1
NrJobsPerBatch <- 200
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
blnL1     <- L1_1000G[,(NrInfoCol + 1):ncol(L1_1000G)] > 0
blnEnough <- rowSums(blnL1, na.rm = T) >= MinNrCarrier
sum(blnEnough)

########################################
#                                      #
#    Run selscan with subset           #
#                                      #
########################################

cat("\n*********  Running selscan  ************\n")

# Get file names, loop over files and do the filtering
SubsetFiles <- list.files(DataFolder, pattern = "NoMulti", full.names = T)

# Set up vectors of start and end indices
Starts <- seq(1, sum(blnEnough), NrJobsPerBatch)
Ends   <- c(Starts[-1] - 1, sum(blnEnough))

# Loop over file names, read file and append to existing
# for (j in 3){
#   cat("\n******   Submitting jobs", Starts[j], "to", Ends[j], "   *********\n")
  for (idxL1 in which(blnEnough)[501:sum(blnEnough)]){
    # Collect chromosome, position, uper and lower bund, and ID for current L1
    Chrom <- paste("chr",  L1_1000G$CHROM[idxL1], sep = "")
    ID    <- L1_1000G$ID[idxL1]
    Pos   <- L1_1000G$POS[idxL1]
    Lower <- Pos - WindowWidth
    Upper <- Pos + WindowWidth
    
    # Create file names
    InFilePattern <- paste(Chrom, "_L1WindowsubsetNoMulti.vcf", sep = "")
    InFile        <- grep(InFilePattern, SubsetFiles, value = T)
    MapFile <- paste(DataFolder, Chrom, "_", ID, "_map", sep = "")
    VcfFile <- paste(DataFolder, Chrom, "_", ID, "_Var.vcf", sep = "")
    OutFile <- gsub("NoMulti", "Selscan", InFile)
    OutFile <- gsub(".vcf", "", OutFile)
    
    
    # Subset vcf file to get variants around current L1
    SubsetCmd <- paste("cut -s -f 1- ", InFile, 
                       " | awk '{if($2 >= ", Lower, "&& $2 <= ", Upper, '){print $0 > "', 
                       VcfFile,'"; print $1, $3, $2, $2 > "', MapFile, '"}}\'', sep = "")
    SelscanCmd    <- paste("./selscan/bin/linux/selscan --ehh",
                           L1_1000G$ID[idxL1], "--vcf",
                           VcfFile, "--map", MapFile, "--out", OutFile)
    Cmds <- c(SubsetCmd, SelscanCmd)
    ScriptName <- paste("L1selscan", Chrom, ID, sep = "_")
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
                             SlurmCommandLines = Cmds, 
                             scriptName = ScriptName) 
  }
#   # Check whether queue is done 
#   QueueFinished <- CheckQueue(MaxNrTrials = MaxNrTrials,
#                               SleepTime   = SleepTime)
#   
# }
