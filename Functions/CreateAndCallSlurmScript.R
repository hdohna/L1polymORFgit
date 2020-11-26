##############################################
#
# General description:
#
#   The following function creates and calls a qsub script

# Input:
#
#     file: the path to save the script file to
#     SlurmHeaderLines: header lines of pbs script
#     SlurmCommandLines: command lines of pbs script
#     scriptName: name to be used to identify run
#     Args: external input arguments

# Output:
#   


# Comments:
#    

##############################################


CreateAndCallSlurmScript <- function(file,
   SlurmHeaderLines = c('#!/bin/bash', 
                       '#SBATCH --account=hb54', 
                       '#SBATCH --time=12:00:00', 
                       '#SBATCH --job-name="My Simple Job."', 
                       '#SBATCH --nodes=1', 
                       '#SBATCH --ntasks=1', 
                       '#SBATCH --cpus-per-task=1', 
                       '#SBATCH --mem=50G' 
   ),
   RunTime = '12:00:00',
   Mem = '50G', 
   SlurmCommandLines, 
   blnSendEmails = F,
   scriptName = 'NoName', 
   Args = "", blnIntern = TRUE,
   blnWait = F,
   blnStopIfNotSubmitted = F){
  
  # Replace name, time and memory in header lines
  SlurmHeaderLines[grep("--job-name=", SlurmHeaderLines)] <- 
    paste('#SBATCH --job-name=', scriptName, sep = '')
  SlurmHeaderLines[grep("--time=", SlurmHeaderLines)] <- 
    paste('#SBATCH --time=', RunTime, sep = '')
  SlurmHeaderLines[grep("#--mem=", SlurmHeaderLines)] <- 
    paste('#SBATCH #SBATCH --mem=', Mem, sep = '')
  
  # Add email options
  if (blnSendEmails){
    SlurmHeaderLines <- c(SlurmHeaderLines,
                          '#SBATCH --mail-user=hb54@aub.edu.lb', 
                          '#SBATCH --mail-type=BEGIN,END,FAIL')
    
  }

  # Save script
  writeLines(c(SlurmHeaderLines, SlurmCommandLines), file)
  
  # Run script
  RunCmd <- paste("sbatch", file, Args)
  RunMessage <- system(command = RunCmd, wait = blnWait, intern = blnIntern)
  
  if (length(grep("Submitted batch job ", RunMessage)) == 1){
    RunID <- strsplit(RunMessage, "Submitted batch job ")[[1]][2]
  } else {
    RunID <- NA
    if(blnStopIfNotSubmitted){
      stop("Job not properly launched")
    }
  }
  list(RunMessage = RunMessage, RunID = RunID)
}