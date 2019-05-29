##############################################
#
# General description:
#
#   The following function checks whether slurm jobs completed

# Input:
#
#     JobIDs:         vector of job IDs (can contain NA)
#     CompletionChar: Character string indicating job completion
#     JobDirectory:   directory where slurm files are

# Output:
#    blnFinished: boolean vector indicating for each job ID whether job has
#                    finished


# Comments:
#    

##############################################


CheckJobCompletion <- function(JobIDs,
                               CompletionChar = "End time:",
                               JobDirectory = "/home/hzudohna/"){
  
  blnNA       <- is.na(JobIDs)
  blnFinished <- rep(NA, length(JobIDs))
  SlurmFiles <- paste(JobDirectory, "slurm-", JobIDs[!blnNA], ".out", sep = "")
  blnFinished[!blnNA] <-  sapply(SlurmFiles, function(x){
    FileLines <- readLines(x)
    length(grep(CompletionChar, FileLines)) > 0
  })
  blnFinished
}