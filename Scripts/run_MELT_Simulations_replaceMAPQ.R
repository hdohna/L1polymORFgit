# The following script sets MAPQ to 0 in unmapped reads in bam files that produced
# a MELT error message and then applies MELT to all bam files.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(Rsamtools)
library(rtracklayer)

# Load previously created objects
load("/labs/dflev/hzudohna/1000Genomes/L1Simulated_MELT.RData")

# Boolean variables for different parts of the workflow
blnRunSim           <- F
blnRunBWA           <- F
blnIdxBam           <- F
blnRunMELT          <- F
blnRunSimAnalysis   <- F
blnRunGroupAnalysis <- T

# Specify run parameters
RunTime      <- '76:00:00'
Mem          <- '200G'
RunTime_MELT <- '2:00:00'
Mem_MELT     <- '50G'

# Path to reference data and for 1000 genomes
#RefPath    <- "D:/L1polymORF/Data/"
RefPath     <- "/labs/dflev/hzudohna/RefSeqData/"
RefFilePath <- "/labs/dflev/hzudohna/RefSeqData/hg19.fa"
Path1000G   <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"
PathGroupMELT   <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT_Group/"
PathScratch   <- "/scratch/users/hzudohna/"

# Load necessary objects
load('/labs/dflev/hzudohna/1000Genomes/GRanges_L1_1000Genomes.RData')
load(paste(RefPath, 'L1RefRanges_hg19.Rdata', sep = ""))

# Get all files with simulated genomes
SimGenomes <- list.files(PathScratch, pattern = "_Haplo1.fa", full.names = T)
SimGenome <- SimGenomes[1]

##################################################################
#                                                                #                             
#         Replace MAPQ values                                    #
#                                                                #                             
##################################################################

# Get IDs that are in the first 50 sample columns but not in 
# SampleIDs
SC50 <- SampleColumns[1:50]
IDsnoMELT <- SC50 [which(!SC50 %in% SampleIDs)]

RunIDs_replace <- NULL

# Loop through the IDs and create new bam files with MAPQ set to zero
for (ID in IDsnoMELT){
  
  # Create file names
  FileBam             <- paste(PathScratch, "hg19_", ID, "_sorted.bam", sep = "")
  FileSam             <- gsub(".bam", ".sam", FileBam)
  FileMAPQreplacedSam <- gsub(".sam", "_MapQreplaced.sam", FileSam)
  FileMAPQreplacedBam <- gsub(".sam", ".bam", FileMAPQreplacedSam)
  FileMAPQRepSortBam <-  gsub(".bam", "_sorted.bam", FileMAPQreplacedBam)
  FileMAPQRepIdxBam <-  gsub(".bam", "_sorted.bam", FileMAPQreplacedBam)
  
  # Sam commands to turn bam file into sam file
  SamBamCmds <- c("module load samtools",
                  paste("samtools view",  FileBam, "-h -o", FileSam))
  
  
  # Awk commands to replace MAPQ by 0 for unmapped reads
  AwkCmds <- c(paste("awk '{if(int(($2 % 8) / 4) == 1){",
                   "print $1, $2, $3, $4, 0, $6, $7, $8, $9, $10, $11;",
                   "} else {",
                   "print $0; } }'",
                   FileSam, ">", FileMAPQreplacedSam),
               "echo 'awk completed'")
  
  # Sam commands to turn sam file into bam file
  SamBamCmds2 <- c("module load samtools",
                   paste("samtools view", FileMAPQreplacedSam, "-h -b -o", 
                         FileMAPQreplacedBam))
  
  # Sam commands to sort and index bam file
  SamIdxCmds <- c("module load samtools",
                  paste("samtools sort", FileMAPQreplacedBam, 
                        "-o", FileMAPQRepSortBam),
                  paste("samtools index", FileMAPQRepSortBam),
                  "echo 'bam file for MELT created'")
  
  # Construct MELT command
  MELTCmds <- c("module load bowtie2",
                paste("java -Xmx2g -jar /labs/dflev/hzudohna/MELTv2.1.5/MELT.jar Single -bamfile",
                      FileMAPQRepSortBam, 
                      "-c 7",
                      "-h /labs/dflev/hzudohna/RefSeqData/hg19.fa", 
                      "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                      "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                      "-b MT/GL000207.1/GL000226.1/GL000229.1/GL000231.1/GL000210.1/GL000239.1/GL000235.1/GL000201.1/GL000247.1/GL000245.1/GL000197.1/GL000203.1/GL000246.1/GL000249.1/GL000196.1/GL000248.1/GL000244.1/GL000238.1/GL000202.1/GL000234.1/GL000232.1/GL000206.1/GL000240.1/GL000236.1/GL000241.1/GL000243.1/GL000242.1/GL000230.1/GL000237.1/GL000233.1/GL000204.1/GL000198.1/GL000208.1/GL000191.1/GL000227.1/GL000228.1/GL000214.1/GL000221.1/GL000209.1/GL000218.1/GL000220.1/GL000213.1/GL000211.1/GL000199.1/GL000217.1/GL000216.1/GL000215.1/GL000205.1/GL000219.1/GL000224.1/GL000223.1/GL000195.1/GL000212.1/GL000222.1/GL000200.1/GL000193.1/GL000194.1/GL000225.1/GL000192.1/NC_007605",
                      "-w", PathGroupMELT))
  
  
  # Create a list of all commands
  AllCmds <- c(SamBamCmds, AwkCmds, SamBamCmds2, SamIdxCmds, MELTCmds)
  
  # Create script name and run script
  ScriptName <- "MAPQReplace_MELTScript"
  RunList <- CreateAndCallSlurmScript(file = ScriptName, 
                                      RunTime = '48:00:00',
                                      Mem     = '200G',
                                      SlurmCommandLines = AllCmds)
  
  cat(RunList$RunMessage, "\n")
  RunIDs_replace <- c(RunIDs_replace, RunList$RunID)
  
}
save(list = "RunIDs_replace", file = "/scratch/users/hzudohna/RunIDsMAPQreplace.RData")

##################################################################
#                                                                #                             
#         Run group MELT on simulated genomes                    #
#                                                                #                             
##################################################################




# Initiailize indicators for finished process
QueueIndivFinished    <- F
QueueGroupFinished    <- F
QueueGenotypeFinished <- F
  
  # Get paths to bam files
  BamFiles <- list.files(PathScratch, pattern = "_sorted.bam", full.names = T)
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
                          "-w", PathGroupMELT))
      
      
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
  # Remove files for failed IndivAnalysis 
  #################
  
  cat("\n*************    Removing files with failes IndivAnalysis    *************\n")
  
  # Check for breaks bam files of size 0 and remove all associated files if
  # there are bam files of size 0
  BrkBamFiles   <- list.files(PathGroupMELT, full.names = T,
                            pattern = "_sorted.LINE1.hum_breaks.sorted.bam")
  BrkBamFiles   <- BrkBamFiles[-grep("bam.bai", BrkBamFiles)]
  BrkBam2Remove <- BrkBamFiles[file.size(BrkBamFiles) == 0]
  for (x in  BrkBam2Remove){
    cat("Removing files associated with", BrkBam2Remove, "\n")
    Split1 <- strsplit(x, "\\/")[[1]]
    RemovePattern <- strsplit(Split1[length(Split1)], "\\.")[[1]][1]
    Files2Remove <- list.files(PathGroupMELT, full.names = T, pattern = RemovePattern)
    file.remove(Files2Remove)
    
  }
  
  
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
                        "-discoverydir", PathGroupMELT,
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-n /labs/dflev/hzudohna/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed",
                         "-w", PathGroupMELT))
    
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
                          "-w", PathGroupMELT,
                          "-p", PathGroupMELT))
      
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
                        "-genotypingdir", PathGroupMELT,
                        "-h", RefFilePath,
                        "-t /labs/dflev/hzudohna/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip",
                        "-w", PathGroupMELT,
                        "-p", PathGroupMELT,
                        "-o", PathGroupMELT))
    
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
  save(list = "JobInfo", file = paste(PathGroupMELT, "L1simulated_MELT_Group_JobInfo.RData",
                                      sep = "")) 



