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


CreateAndCallqsubScript_scg4 <- function(file,
   qsubHeaderLines = c('#! /bin/sh', '#', '#$ -N TEST', '#', '#$ -cwd', '#', 
                       '#$ -l h_rt=72:00:00', '#', '#$ -j y', '#',
                       '#$ -P large_mem', '#',
                       '#$ -S /bin/bash', '#', ''), 
   qsubCommandLines, scriptName = 'NoName', Args = ""){
  
  # Replace name in header lines
  qsubHeaderLines[grep("-N", qsubHeaderLines)] <- paste('#$ -N', scriptName)
  
  # Save script
  writeLines(c(qsubHeaderLines, qsubCommandLines), file)
  
  # Run script
  RunCmd <- paste("qsub", file, Args)
  system(RunCmd, wait = F)
  
}