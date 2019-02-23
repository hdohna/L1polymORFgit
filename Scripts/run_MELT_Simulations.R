# The following script creates objects and commands for calling L1 deletions in
# 1000 genome data using MELT.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Boolean variables for different parts of the workflow
blnRunSim  <- F
blnRunBWA  <- F
blnIdxBam  <- F
blnRunMELT <- T
blnRunSimAnalysis   <- F
blnRunGroupAnalysis <- T

# Specify run parameters
RunTime <- '12:00:00'
Mem     <- '100G'
RunTime_MELT <- '6:00:00'
Mem_MELT     <- '50G'

# Path to reference data and for 1000 genomes
#RefPath    <- "D:/L1polymORF/Data/"
RefPath     <- "/labs/dflev/hzudohna/RefSeqData/"
RefFilePath <- "/labs/dflev/hzudohna/RefSeqData/hg19.fa"
Path1000G   <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"

# Load necessary objects
load(paste(Path1000G, 'GRanges_L1_1000Genomes.RData', sep = ""))
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Get all files with simulated genomes
SimGenomes <- list.files(Path1000G, pattern = "_Haplo1.fa", full.names = T)
SimGenome <- SimGenomes[1]

##################################################################
#                                                                #                             
#            Simulate and align reads for genome,                #
#             run single MELT on them                            #
#                                                                #                             
##################################################################

if (blnRunSimAnalysis){
  for (SimGenome in SimGenomes){
    
    # Get ID
    Split1 <- strsplit(SimGenome, "_")[[1]]
    ID     <- strsplit(Split1[length(Split1)], "\\.")[[1]][1]
    
    # Create folder for MELT
    NewDir <- paste(Path1000G, Split1[4], sep = "")
    
    # Generate output file names for simulation 
    SimGenome2  <- gsub("1.fa", "2.fa", SimGenome)
    
    # Generate output file names for simulation 
    Fq1     <- gsub("_Haplo1.fa", "_1.fq", SimGenome)
    Fq11    <- gsub("_Haplo1.fa", "_11.fq", SimGenome)
    Fq12    <- gsub("_Haplo1.fa", "_12.fq", SimGenome)
    Fq2     <- gsub("_Haplo1.fa", "_2.fq", SimGenome)
    Fq21    <- gsub("_Haplo1.fa", "_21.fq", SimGenome)
    Fq22    <- gsub("_Haplo1.fa", "_22.fq", SimGenome)
    Sai1   <- gsub("_Haplo1.fa", "_1.sai", SimGenome)
    Sai2   <- gsub("_Haplo1.fa", "_2.sai", SimGenome)
    OutSim1 <- gsub("_Haplo1.fa", "_wgsim_out1", SimGenome)
    OutSim2 <- gsub("_Haplo1.fa", "_wgsim_out2", SimGenome)
    OutSam  <- gsub("_Haplo1.fa", ".sam", SimGenome)
    OutBam  <- gsub("_Haplo1.fa", ".bam", SimGenome)
    OutBamSorted <- gsub("_Haplo1.fa", "_sorted.bam", SimGenome)
    OutVcf <- gsub("_Haplo1.fa", "_sorted.vcf", SimGenome)
    
    # Construct command to simulate genome
    wgsimCmds <- c(paste("/home/hzudohna/wgsim/wgsim -h -e 0.004 -d 379.5 -s 21.6 -1 100 -2 100 -N 70000000",
                         SimGenome, Fq11, Fq12, ">", OutSim1),
                   paste("/home/hzudohna/wgsim/wgsim -h -e 0.004 -d 379.5 -s 21.6 -1 100 -2 100 -N 70000000",
                         SimGenome2, Fq21, Fq22, ">", OutSim2),
                   paste("cat", Fq11, Fq21, ">", Fq1),
                   paste("cat", Fq12, Fq22, ">", Fq2),
                   paste("rm", Fq11, Fq21, Fq12, Fq22))
    
    
    # Bwa commands
    BWACmds <- c("module load bwa",
                 paste("bwa aln", RefFilePath, Fq1, 
                       paste(">", Sai1)),
                 paste("bwa aln", RefFilePath, Fq2, 
                       paste(">", Sai2)),
                 paste("bwa sampe -a 618", RefFilePath,
                       Sai1, Sai2, Fq1, Fq2, ">", OutSam))
    
    # Sam commands to turn sam file into bam file
    SamBamCmds <- c("module load samtools",
                    paste("samtools view", OutSam, "-b -h -o", OutBam),
                    paste("rm", OutSam))
    
    # Sam commands to turn sam file into bam file
    SamIdxCmds <- c("module load samtools",
                    paste("samtools sort", OutBam, "-o", OutBamSorted),
                    paste("samtools index", OutBamSorted))
    
    # Construct command to make directory
    MkDirCmd <- NULL
    if (!dir.exists(NewDir)){
      paste("mkdir", NewDir)
    }
    
    # Construct MELT command
    MELTCmds <- c(MkDirCmd,
                  "module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
                        OutBamSorted, 
                        "-c 7",
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                        "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                        "-w", NewDir))
    
    
    # Create a list of all commands
    CmdList <- list(wgsimCmds, c(BWACmds, SamBamCmds), SamIdxCmds, MELTCmds)
    AllCmds <- unlist(CmdList[c(blnRunSim, blnRunBWA, blnIdxBam, blnRunMELT)])
    
    # Create script name and run script
    ScriptName <- paste("Sim_Bwa_MELTScript", ID, sep = "_")
    CreateAndCallSlurmScript(file = ScriptName, 
                             RunTime = RunTime,
                             Mem = Mem,
                             SlurmCommandLines = AllCmds)
    
  }
  
}

##################################################################
#                                                                #                             
#         Run group MELT on simulated genomes                    #
#                                                                #                             
##################################################################

if(blnRunGroupAnalysis){
  
  # Get paths to bam files
  BamFiles <- list.files(Path1000G, pattern = "_sorted.bam", full.names = T)
  BamFiles <- BamFiles[-grep("_sorted.bam.", BamFiles)]
  
  # Perform variant discovery for each individual bam file
  cat("\n*************    Performing MELT for individual bam files     *************\n")
  for (BamFile in BamFiles){
    
    # Get ID from bam file
    Split1 <- strsplit(BamFile, "_")[[1]]
    ID     <- strsplit(Split1[length(Split1) - 1], "\\.")[[1]][1]
    cat("Analyzing", ID, "\n")
    
    # Create commands to run MELT
    MELTCmds <- c("module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar IndivAnalysis",
                        " -bamfile", BamFile, 
                        "-c 7",
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                        "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                        "-w", Path1000G))
    
    
    # Create script name and run script
    ScriptName <- paste("GroupInd_MELTScript", ID, sep = "_")
    CreateAndCallSlurmScript(file = ScriptName, 
                             RunTime = RunTime,
                             Mem = Mem,
                             SlurmCommandLines = MELTCmds)
    
  }
  
  # Check whether for bam files queue is finished, once it is finished, run group analysis
  QueueBamFinished <- CheckQueue(MaxNrTrials = 30, SleepTime   = 30)
  if (QueueBamFinished){
    
    cat("\n*************   Summarizing MELT for individual bam files     *************\n")
    # Create directory for discovery file
    DiscoveryDir <- paste(Path1000G, "L1DiscoveryGroup", sep = "") 
    
    # Construct command to make directory
    MkDiscoveryDirCmd <- NULL
    if (!dir.exists(DiscoveryDir)){
      MkDiscoveryDirCmd <- paste("mkdir", DiscoveryDir)
    }

    # Create commands to run MELT
    MELTCmds <- c(MkDiscoveryDirCmd,
                  "module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar GroupAnalysis",
                        "-discoverydir", DiscoveryDir,
                        "-c 7",
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                        "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                        "-w", Path1000G))
    
    
    # Create script name and run script
    ScriptName <- paste("GroupInd_MELTScript", ID, sep = "_")
    CreateAndCallSlurmScript(file = ScriptName, 
                             RunTime = RunTime,
                             Mem = Mem,
                             SlurmCommandLines = MELTCmds)
    
    # Check whether queue of group analysis is finished, once it is finished, run group analysis
    QueueGroupFinished <- CheckQueue(MaxNrTrials = 30, SleepTime   = 30)
    if (QueueGroupFinished){
      
      cat("\n*************   Getting genotypes for individual bam files     *************\n")
      
      # Perform variant discovery for each individual bam file
      for (BamFile in BamFiles){
        
        # Get ID from bam file
        Split1 <- strsplit(BamFile, "_")[[1]]
        ID     <- strsplit(Split1[length(Split1)], "\\.")[[1]][1]
        
        # Create commands to run MELT
        MELTCmds <- c("module load bowtie2",
                      paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Genotype",
                            " -bamfile", BamFile, 
                            "-c 7",
                            "-h", RefFilePath,
                            "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                            "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                            "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                            "-w", DiscoveryDir))
        
        
        # Create script name and run script
        ScriptName <- paste("GroupInd_MELTScript", ID, sep = "_")
        CreateAndCallSlurmScript(file = ScriptName, 
                                 RunTime = RunTime,
                                 Mem = Mem,
                                 SlurmCommandLines = MELTCmds)
        
      }
    }
    
    
    
  }
}


