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
Path1000G <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Get all files with simulated genomes
SimGenomes <- list.files(Path1000G, pattern = ".fa", full.names = T)

########################################################################################
# Simulate and align reads for enome
SimGenome <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT//hg19_HG00096.fa"

for (SimGenome in SimGenomes[-1]){
  
  # Get ID
  Split1 <- strsplit(SimGenome, "_")[[1]]
  ID     <- strsplit(Split1[length(Split1)], "\\.")[[1]][1]
  
  # Generate output file names for simulation 
  Fq1    <- gsub(".fa", "_1.fq", SimGenome)
  Fq2    <- gsub(".fa", "_2.fq", SimGenome)
  Sai1   <- gsub(".fa", "_1.sai", SimGenome)
  Sai2   <- gsub(".fa", "_2.sai", SimGenome)
  OutSim <- gsub(".fa", "_wgsim_out", SimGenome)
  OutSam <- gsub(".fa", ".sam", SimGenome)
  
  # Construct MELT command
  wgsimCmds <- paste("/home/hzudohna/wgsim/wgsim -e 0.004 -d 379.5 -s 21.6 -1 100 -2 100 -N 145000000",
                 SimGenome, Fq1, Fq2, ">", OutSim)
  
  # Bwa commands
  BWACmds <- c("module load bwa",
               paste("bwa aln /labs/dflev/hzudohna/RefSeqData/hg19.fa", Fq1, 
                     paste(">", Sai1)),
               paste("bwa aln /labs/dflev/hzudohna/RefSeqData/hg19.fa", Fq2, 
                     paste(">", Sai2)),
               paste("bwa sampe -a 618 /reference/RefGenomes/H_sapiens/hg19/hg19.fa",
                     Sai1, Sai2, Fq1, Fq2, ">", OutSam))

  # Create script name and run script
  ScriptName <- paste("Sim_Bwa_Script", ID, sep = "_")
  CreateAndCallSlurmScript(file = ScriptName, 
                           SlurmCommandLines = c(wgsimCmds, 
                                                 "echo 'simulation completed'", 
                                                 BWACmds,
                                                 "echo 'alignment completed'"))
  
}

