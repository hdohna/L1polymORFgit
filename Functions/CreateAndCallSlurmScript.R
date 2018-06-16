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
                       '#SBATCH --account=dflev', 
                       '#SBATCH --time=1:00:00', 
                       '#SBATCH --job-name="My Simple Job."', 
                       '#SBATCH --nodes=1', 
                       '#SBATCH --ntasks=1', 
                       '#SBATCH --cpus-per-task=1', 
                       '#SBATCH --mem=4G', 
                       '#SBATCH --mail-user=hb54@aub.edu.lb', 
                       '#SBATCH --mail-type=BEGIN,END,FAIL'
   ), 
   SlurmCommandLines, 
   scriptName = 'NoName', 
   Args = "",
   blnWait = F){
  
  # Replace name in header lines
  SlurmHeaderLines[grep("--job-name=", SlurmHeaderLines)] <- 
    paste('#SBATCH --job-name=', scriptName)
  
  # Save script
  writeLines(c(SlurmHeaderLines, SlurmCommandLines), file)
  
  # Run script
  RunCmd <- paste("sbatch", file, Args)
  system(RunCmd, wait = blnWait)
  
}