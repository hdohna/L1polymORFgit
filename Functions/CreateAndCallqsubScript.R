##############################################
#
# General description:
#
#   The following function creates and calls a qsub script

# Input:
#
#     file: the path to save the script file to
#     qsubHeaderLines: header lines of pbs script
#     qsubCommandLines: command lines of pbs script
#     scriptName: name to be used to identify run
#     Args: external input arguments

# Output:
#   


# Comments:
#    

##############################################


CreateAndCallqsubScript <- function(file,
   qsubHeaderLines = c('#!/bin/bash', 
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
   qsubCommandLines, 
   scriptName = 'NoName', 
   Args = "",
   blnWait = F){
  
  # Replace name in header lines
  qsubHeaderLines[grep("--job-name=", qsubHeaderLines)] <- 
    paste('#SBATCH --job-name=', scriptName)
  
  # Save script
  writeLines(c(qsubHeaderLines, qsubCommandLines), file)
  
  # Run script
  RunCmd <- paste("sbatch", file, Args)
  system(RunCmd, wait = blnWait)
  
}