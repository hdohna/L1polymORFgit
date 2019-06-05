# The following script creates objects and commands for calling L1 deletions in
# 1000 genome data using MELT.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Boolean variables for different parts of the workflow
blnRunSim           <- T
blnRunBWA           <- T
blnIdxBam           <- T
blnRunMELT          <- T
blnRunSimAnalysis   <- T
blnRunGroupAnalysis <- T

blnRunSim_Var           <- F
blnRunBWA_Var           <- F
blnIdxBam_Var           <- F
blnRunMELT_Var          <- F
blnRunSimAnalysis_Var   <- F
blnRunGroupAnalysis_Var <- F

# Specify run parameters
RunTime      <- '24:00:00'
Mem          <- '200G'
RunTime_MELT <- '2:00:00'
Mem_MELT     <- '50G'

# Path to reference data and for 1000 genomes
#RefPath    <- "D:/L1polymORF/Data/"
RefPath     <- "/labs/dflev/hzudohna/RefSeqData/"
RefFilePath <- "/labs/dflev/hzudohna/RefSeqData/hg19.fa"
Path1000G   <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"
PathScratch   <- "/scratch/users/hzudohna/"

# Load necessary objects
load('/labs/dflev/hzudohna/1000Genomes/GRanges_L1_1000Genomes.RData')
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
  cat("***********  Aligning simulated genomes    ******* \n")
  for (SimGenome in SimGenomes){
    
    # Get ID
    Split1 <- strsplit(SimGenome, "_")[[1]]
    ID     <- strsplit(Split1[length(Split1) - 1], "\\.")[[1]][1]
    
    # Create folder for MELT
    NewDir <- paste(Path1000G, Split1[4], sep = "")
    
    # Generate output file names for simulation 
    SimGenome2  <- gsub("1.fa", "2.fa", SimGenome)
    
    # Generate a temporary simulation path
    SimGenome_tmp <- gsub(Path1000G, PathScratch, SimGenome)
    
    # Generate output file names for simulation 
    Fq1     <- gsub("_Haplo1.fa", "_1.fq", SimGenome_tmp)
    Fq11    <- gsub("_Haplo1.fa", "_11.fq", SimGenome_tmp)
    Fq12    <- gsub("_Haplo1.fa", "_12.fq", SimGenome_tmp)
    Fq2     <- gsub("_Haplo1.fa", "_2.fq", SimGenome_tmp)
    Fq21    <- gsub("_Haplo1.fa", "_21.fq", SimGenome_tmp)
    Fq22    <- gsub("_Haplo1.fa", "_22.fq", SimGenome_tmp)
    Sai1   <- gsub("_Haplo1.fa", "_1.sai", SimGenome_tmp)
    Sai2   <- gsub("_Haplo1.fa", "_2.sai", SimGenome_tmp)
    OutSim1 <- gsub("_Haplo1.fa", "_wgsim_out1", SimGenome_tmp)
    OutSim2 <- gsub("_Haplo1.fa", "_wgsim_out2", SimGenome_tmp)
    OutSam  <- gsub("_Haplo1.fa", ".sam", SimGenome_tmp)
    OutBam  <- gsub("_Haplo1.fa", ".bam", SimGenome_tmp)
    OutBamSorted <- gsub("_Haplo1.fa", "_sorted.bam", SimGenome_tmp)
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
                    paste("rm", OutSam),
                    paste("rm", Fq1),
                    paste("rm", Fq2))
    
    # Sam commands to sort and index bam file
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
#         Simulate and align reads for genomes for               #
#          variable differences from consensus L1,               #
#             run single MELT on them                            #
#                                                                #                             
##################################################################

# Get all files with simulated genomes
SimGenomes <- list.files(Path1000G, pattern = "_Haplo1.fa", full.names = T)
SimGenomes <- grep("hg19Var_", SimGenomes, value = T)
SimGenome <- SimGenomes[1]

if (blnRunSimAnalysis_Var){
  cat("***********  Running simulations for L1 with variable differences from consensus  *******\n")
  for (SimGenome in SimGenomes[1:20]){
    
    # Get ID
    Split1 <- strsplit(SimGenome, "_")[[1]]
    ID     <- strsplit(Split1[length(Split1) - 1], "\\.")[[1]][1]
    
    # Create folder for MELT
    NewDir <- paste(Path1000G, Split1[4], "_Var", sep = "")
    
    # Generate output file names for simulation 
    SimGenome2  <- gsub("1.fa", "2.fa", SimGenome)
    
    # Generate a temporary simulation path
    SimGenome_tmp <- gsub(Path1000G, PathScratch, SimGenome)

    # Generate output file names for simulation 
    Fq1     <- gsub("_Haplo1.fa", "_1.fq", SimGenome)
    Fq11    <- gsub("_Haplo1.fa", "_11.fq", SimGenome)
    Fq12    <- gsub("_Haplo1.fa", "_12.fq", SimGenome)
    Fq2     <- gsub("_Haplo1.fa", "_2.fq", SimGenome)
    Fq21    <- gsub("_Haplo1.fa", "_21.fq", SimGenome)
    Fq22    <- gsub("_Haplo1.fa", "_22.fq", SimGenome)
    Sai1    <- gsub("_Haplo1.fa", "_1.sai", SimGenome)
    Sai2    <- gsub("_Haplo1.fa", "_2.sai", SimGenome)
    OutSim1 <- gsub("_Haplo1.fa", "_wgsim_out1", SimGenome_tmp)
    OutSim2 <- gsub("_Haplo1.fa", "_wgsim_out2", SimGenome_tmp)
    OutSam  <- gsub("_Haplo1.fa", ".sam", SimGenome_tmp)
    OutBam  <- gsub("_Haplo1.fa", ".bam", SimGenome_tmp)
    OutBamSorted <- gsub("_Haplo1.fa", "_sorted.bam", SimGenome_tmp)
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
                    paste("rm", OutSam),
                    paste("rm", Fq1),
                    paste("rm", Fq2))

    # Sam commands to sort and index bam file
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
    AllCmds <- unlist(CmdList[c(blnRunSim_Var, blnRunBWA_Var, blnIdxBam_Var, blnRunMELT_Var)])
    
    # Create script name and run script
    ScriptName <- paste("Sim_Bwa_MELTScript", ID, sep = "_")
    CreateAndCallSlurmScript(file = ScriptName, 
                             scriptName = paste(ID, "VarSim", sep = "_"),
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
  
  # Initiailize indicators for finished process
  QueueIndivFinished    <- F
  QueueGroupFinished    <- F
  QueueGenotypeFinished <- F
  
  # Get paths to bam files
  BamFiles <- list.files(Path1000G, pattern = "_sorted.bam", full.names = T)
  BamFiles <- BamFiles[-grep("_sorted.bam.", BamFiles)]
  
  #################
  # Preprocess individual bam file
  #################

  cat("\n*************    Preprocessing individual bam files     *************\n")
  RunID <- NA
  RunIDs_Preprocess <- NULL
  for (BamFile in BamFiles){
    
    # Get ID from bam file
    Split1 <- strsplit(BamFile, "_")[[1]]
    ID     <- strsplit(Split1[length(Split1) - 1], "\\.")[[1]][1]
    cat("Analyzing", ID, "\n")
    
    # Create commands to run MELT
    MELTCmds <- c("module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Preprocess",
                        "-bamfile", BamFile,
                        "-h", RefFilePath
                  ))
    
    # Create script name and run script
    ScriptName <- paste("GroupInd_Preprocess", ID, sep = "_")
    RunID <- CreateAndCallSlurmScript(file = ScriptName, 
                               RunTime = RunTime_MELT,
                               Mem = Mem_MELT,
                               SlurmCommandLines = MELTCmds)$RunID
    RunIDs_Preprocess <- c(RunIDs_Preprocess, RunID)
  }
  
  # Check whether for bam files queue is finished, once it is finished, run group analysis
  cat("\n\n")
  Sys.sleep(10)
  QueuePreprocessFinished <- CheckQueue(MaxNrTrials = 500, SleepTime = 60,
                                        JobIDs = RunIDs_Preprocess)
  blnFinished <- CheckJobCompletion(RunIDs_Preprocess)
  JobInfo     <- data.frame(JobID = RunIDs_Preprocess, JobType = "MELT_Preprocess", 
                        InputFile = BamFiles, blnFinished  = blnFinished)
    
  #################
  # Perform IndivAnalysis (variant discovery for each individual bam file)
  #################
  
  RunIDs_IndivAnalysis  <- NULL
  blnFinished <- NULL
  RunID <- NA
  if(QueuePreprocessFinished){
    cat("\n*************    Performing MELT IndivAnalysis     *************\n")
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
                          "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                          "-w", Path1000G))
      
      
      # Create script name and run script
      ScriptName <- paste("GroupInd_MELTScript", ID, sep = "_")
      RunID <- CreateAndCallSlurmScript(file = ScriptName, 
                               RunTime = RunTime_MELT,
                               Mem = Mem_MELT,
                               SlurmCommandLines = MELTCmds)$RunID
      RunIDs_IndivAnalysis <- c(RunIDs_IndivAnalysis, RunID)
      
    } # End of loop over bam files
    
    # Check whether for bam files queue is finished, once it is finished, run group analysis
    cat("\n\n")
    Sys.sleep(10)
    QueueIndivFinished <- CheckQueue(MaxNrTrials = 50, SleepTime = 30,
                                     JobIDs = RunIDs_IndivAnalysis)
    blnFinished <- CheckJobCompletion(RunIDs_IndivAnalysis)
    NewJobInfo <- data.frame(JobID = RunIDs_IndivAnalysis, JobType = "MELT_IndivAnalysis", 
                          InputFile = BamFiles, blnFinished  = blnFinished)
    JobInfo <- rbind(JobInfo, NewJobInfo)
    
    
  } # End of IndivAnalysis
  
  
  #################
  # Perform GroupAnalysis
  #################

  RunID <- NA
  RunIDs_GroupAnalysis <- NULL
  if (QueueIndivFinished){
    
    cat("\n*************   Performing MELT GroupAnalysis      *************\n")
    # Create commands to run MELT
    MELTCmds <- c("module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar GroupAnalysis",
                        "-discoverydir", Path1000G,
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                         "-w", Path1000G))
    
    # Create script name and run script
    ScriptName <- paste("GroupAnalysis_MELTScript")
    RunID <- CreateAndCallSlurmScript(file = ScriptName, 
                             RunTime = RunTime_MELT,
                             Mem = Mem_MELT,
                             SlurmCommandLines = MELTCmds)$RunID
    RunIDs_GroupAnalysis <- c(RunIDs_GroupAnalysis, RunID)
    
    # Check whether queue of group analysis is finished, once it is finished, run group analysis
    Sys.sleep(10)
    QueueGroupFinished <- CheckQueue(MaxNrTrials = 50, SleepTime = 30,
                                     JobIDs = RunIDs_GroupAnalysis)
    blnFinished <- CheckJobCompletion(RunIDs_GroupAnalysis)
    NewJobInfo <- data.frame(JobID = RunIDs_GroupAnalysis, JobType = "MELT_GroupAnalysis", 
                             InputFile = NA, blnFinished  = blnFinished)
    JobInfo <- rbind(JobInfo, NewJobInfo)
    
    
  } # End of GroupAnalysis
  
  #################
  # Perform Genotype
  #################
  
  RunID      <- NA
  RunIDs_Genotype <- NULL
  if (QueueGroupFinished){
    
    cat("\n*************   Getting genotypes for individual bam files     *************\n")
    for (BamFile in BamFiles){
      
      # Get ID from bam file
      Split1 <- strsplit(BamFile, "_")[[1]]
      ID     <- strsplit(Split1[length(Split1) - 1], "\\.")[[1]][1]
      cat("Analyzing", ID, "\n")
      
      # Create commands to run MELT
      MELTCmds <- c("module load bowtie2",
                    paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Genotype",
                          " -bamfile", BamFile, 
                          "-h", RefFilePath,
                          "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                          "-w", Path1000G,
                          "-p", Path1000G))
      
      # Create script name and run script
      ScriptName <- paste("Genotype_MELTScript", ID, sep = "_")
      RunID <- CreateAndCallSlurmScript(file = ScriptName, 
                               RunTime = RunTime_MELT,
                               Mem = Mem_MELT,
                               SlurmCommandLines = MELTCmds)$RunID
      RunIDs_Genotype <- c(RunIDs_Genotype, RunID)
      
    } # End of loop over bam files to genotype
    
    # check whether genotyping loop has finished
    cat("\n\n")
    Sys.sleep(10)
    QueueGenotypeFinished <- CheckQueue(MaxNrTrials = 30, SleepTime = 30,
                                        JobIDs = RunIDs_Genotype)
    blnFinished <- CheckJobCompletion(RunIDs_Genotype)
    NewJobInfo <- data.frame(JobID = RunIDs_Genotype, JobType = "MELT_Genotype", 
                             InputFile = BamFiles, blnFinished  = blnFinished)
    JobInfo <- rbind(JobInfo, NewJobInfo)
    
  } # End of MELT Genotype
  
  #################
  # Perform MakeVCF
  #################

  RunID     <- NA
  RunIDs_MakeVCF <- NULL
  if (QueueGenotypeFinished){
    
    cat("\n*************   Run MELT MakeVCF     *************\n")
    # Create commands to run MELT
    MELTCmds <- c("module load bowtie2",
                  paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar MakeVCF",
                        "-genotypingdir", Path1000G,
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-w", Path1000G,
                        "-p", Path1000G,
                        "-o", Path1000G))
    
    # Create script name and run script
    ScriptName <- paste("Vcf_MELTScript")
    RunID <- CreateAndCallSlurmScript(file = ScriptName, 
                             RunTime = RunTime_MELT,
                             Mem = Mem_MELT,
                             SlurmCommandLines = MELTCmds)$RunID
    RunIDs_MakeVCF <- c(RunIDs_MakeVCF, RunID)
    cat("MakeVcf is job", RunIDs_MakeVCF, "\n")
    QueueMakeVCFFinished <- CheckQueue(MaxNrTrials = 30, SleepTime = 30,
                                        JobIDs = RunIDs_MakeVCF)
    blnFinished <- CheckJobCompletion(RunIDs_MakeVCF)
    NewJobInfo <- data.frame(JobID = RunIDs_MakeVCF, JobType = "MELT_MakeVCF", 
                             InputFile = NA, blnFinished  = blnFinished)
    JobInfo <- rbind(JobInfo, NewJobInfo)
    
  } # End of MELT MELT_MakeVCF
  
  # Code that can be used to check status of individual jobs
  # JobIDs =   JobInfo$JobID[JobInfo$JobType == "MELT_MakeVCF"]
  # SlurmFiles <- paste(JobDirectory, "slurm-", JobIDs, ".out", sep = "")
  # lapply(SlurmFiles, function(x) readLines(x))
  
  # Save job info
  save(list = "JobInfo", file = paste(Path1000G, "L1simulated_MELT_Group_JobInfo.RData",
                                      sep = "")) 
}


