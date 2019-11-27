# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

RunIDs_sample <- NULL
for (i in 1:100){
  ScriptCmd <- c("module load r",
                 paste("Rscript /home/hzudohna/L1polymORFgit/Scripts/EstimateL1SelectionPars_MELT_Group_SingleSample.R",
                       i))
  ScriptName <- paste("L1SelectionSample", i, sep = "_")
  RunList <- CreateAndCallSlurmScript(file = ScriptName, 
                                      RunTime = '10:00:00',
                                      Mem     = '50G',
                                      SlurmCommandLines = ScriptCmd)
  
  cat(RunList$RunMessage, "\n")
  RunIDs_sample <- c(RunIDs_sample, RunList$RunID)
  
}
save(list = "RunIDs_sample", 
     file = "/labs/dflev/hzudohna/1000Genomes/L1Selection_Sample_runIDs.RData")