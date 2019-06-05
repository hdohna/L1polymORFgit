# The following script creates a reference genome to test LINE-1 detection
load('/labs/dflev/hzudohna/RefSeqData/GRanges_L1_1000Genomes.RData')

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Set parameters
Samples2Use <- SampleColumns[1:50]


##################################################################
#                                                                #                             
#              Create genomes with variable                      #
#                 L1 insertion length                            #
#                                                                #                             
##################################################################

# Loop through the first sample columns and generate reference genomes
# with the same insertions as the current
ScriptPathGeneric <- "/scratch/users/hzudohna/CreateTestGenome_Generic.R"
RunIDs_CreateGenomes <- NULL

for (x in Samples2Use){
  cat("******   Simulating genome", x, "   **********\n")
  ScriptLines <- readLines('/home/hzudohna/L1polymORFgit/Scripts/CreateTestGenome_Generic.R')
  ScriptLines <- gsub("GenericID", x, ScriptLines)
  OutPath     <- gsub("Generic", x, ScriptPathGeneric)
  writeLines(ScriptLines, con = OutPath)
  Cmd <- paste("sbatch sb_R", OutPath)
  RunMessage <- system(Cmd, intern = T)
  cat(RunMessage, "\n")
  RunID <- strsplit(RunMessage, "Submitted batch job ")[[1]][2]
  RunIDs_CreateGenomes <- c(RunIDs_CreateGenomes, RunID)
  
}


##################################################################
#                                                                #                             
#              Check that genomes have the                       #
#              specified L1 sequences                            #
#                                                                #                             
##################################################################

# Check whether creating genomes has finished and if so, check the L1
# sequences
cat("**********   Checking whether creation of genomes has finished     *********\n")
QueueCreateGenomesFinished <- CheckQueue(MaxNrTrials = 100, SleepTime = 30,
                                      JobIDs = RunIDs_CreateGenomes)
ScriptPathGeneric <- "/scratch/users/hzudohna/CheckTestGenome_Generic.R"
RunIDs_CheckGenomes <- NULL

# Check whether L1 sequences are as expected when the queue has finished
if (QueueCreateGenomesFinished){
  for (x in Samples2Use){
    cat("******   Checking genome", x, "   **********\n")
    ScriptLines <- readLines('/home/hzudohna/L1polymORFgit/Scripts/CheckTestGenome_Generic.R')
    ScriptLines <- gsub("GenericID", x, ScriptLines)
    OutPath     <- gsub("Generic", x, ScriptPathGeneric)
    writeLines(ScriptLines, con = OutPath)
    Cmd <- paste("sbatch sb_R", OutPath)
    RunMessage <- system(Cmd, intern = T)
    cat(RunMessage, "\n")
    RunID <- strsplit(RunMessage, "Submitted batch job ")[[1]][2]
    RunIDs_CheckGenomes <- c(RunIDs_CheckGenomes, RunID)
  }
}

# Check whether checking has finished and if so, return the status
cat("**********   Checking whether checking of genomes has finished     *********\n")
QueueCheckGenomesFinished <- CheckQueue(MaxNrTrials = 100, SleepTime = 30,
                                         JobIDs = RunIDs_CheckGenomes)

if (QueueCheckGenomesFinished){
  JobDirectory <- getwd()
  SlurmFiles <- paste(JobDirectory, "/slurm-", RunIDs_CheckGenomes, ".out", sep = "")
  SlurmFiles <- SlurmFiles[file.exists(SlurmFiles)]
  L1Status <- sapply(SlurmFiles, function(x){
    FileLines <- readLines(x)
    c(grep("Total match:", FileLines, value = T),
      grep("Total mismatch:", FileLines, value = T))
  })
  
}

# Create an overview table and write it out
OverviewTable <- data.frame(SampleID = Samples2Use,
                            GCreationID = RunIDs_CreateGenomes,
                            GCheckID    = RunIDs_CheckGenomes,
                            L1Match = L1Status[1, ],
                            L1MisMatch = L1Status[2, ])
outPath <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/SimGenomesOverview.csv"
cat("Writing overview table to", outPath, ".....")
write.csv(OverviewTable, outPath)
cat("done!\n")

##################################################################
#                                                                #                             
#              Create genomes with variable                      #
#               difference from consensus                        #
#                                                                #                             
##################################################################

# # Loop through the first sample columns and generate reference genomes
# # with the same insertions as the current
# for (x in SampleColumns[1:20]){
#   
#   cat("******   Simulating genome", x, "   **********\n")
#   
#   # Path for genome fasta file path
#   FastaPath1 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19Var_",
#                       x, "_Haplo1.fa", sep = "")
#   FastaPath2 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19Var_",
#                       x, "_Haplo2.fa", sep = "")
#   if (file.exists(FastaPath1)){
#     file.remove(FastaPath1)
#   }
#   if (file.exists(FastaPath2)){
#     file.remove(FastaPath2)
#   }
#   
#   # Open connections for fasta file to write chromosomes
#   GenCon1 <- file(FastaPath1, open = "a")
#   GenCon2 <- file(FastaPath2, open = "a")
#   
#   # Index of L1 in current genome and all chromosomes with L1
#   idxL1   <- which(L1_1000G[,x] > 0)
#   ChrNrs <- L1_1000G$CHROM[idxL1]
#   Chroms <- paste("chr", ChrNrs, sep = "")
#   UniqueChroms <- unique(Chroms)
#   
#   # Loop over chromosomes and generate insertions
#   for (Chr in AllChrs){
#     cat("Processing", Chr, "\n")
#     CurrentChrom  <- paste('BSgenome.Hsapiens.UCSC.hg19[["', Chr, '"]]', sep = "")
#     NewDNASt_txt1 <- CurrentChrom
#     NewDNASt_txt2 <- CurrentChrom
#     if(Chr %in% UniqueChroms){
#       idxChr     <- idxL1[Chr == Chroms]
#       Count      <- L1_1000G[idxChr,x]
#       bln2       <- Count == 2
#       
#       # Create indices for both homologous chromosome
#       idx1       <- idxChr[bln2]
#       idx2       <- idxChr[bln2]
#       HaploSample <- sample(c(T, F), sum(!bln2), replace = T)
#       idx1 <- c(idx1, idxChr[which(!bln2)[HaploSample]])
#       idx2 <- c(idx2, idxChr[which(!bln2)[!HaploSample]])
#       
#       # Create insertion patterns for both homologous chromosome
#       if (length(idx1) > 0) NewDNASt_txt1 <- CreateInsertTxtVar(idx1)
#       if (length(idx2) > 0) NewDNASt_txt2 <- CreateInsertTxtVar(idx2)
#     }
#     NewDNASt1      <- eval(parse(text = NewDNASt_txt1))
#     NewDNASt2      <- eval(parse(text = NewDNASt_txt2))
#     NewDNASt_char1 <- as.character(NewDNASt1)
#     NewDNASt_char2 <- as.character(NewDNASt2)
#     ChromName <- paste(">", substr(Chr, 4, nchar(Chr)), sep = "")
#     writeLines(text = c(ChromName, NewDNASt_char1), con = GenCon1)
#     writeLines(text = c(ChromName, NewDNASt_char2), con = GenCon2)
#   }
#   close(GenCon1)  
#   close(GenCon2)  
#   
# }
# 
# # Save settings
# save.image(file = "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/L1SimSettings.RData")