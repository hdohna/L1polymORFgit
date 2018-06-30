##############################################
#
# General description:
#
#   The following checks a slurm queue repeatedly and
#   returns a boolean whether the queue is finished

# Input:
#
#     MaxNrTrials: integer for maximum number of trials
#     SleepTime: how long to wait until next check

# Output:
#    QueueFinished: logical indicating whether queue ha finished


# Comments:
#    

##############################################


CheckQueue <- function(MaxNrTrials = 10,
                       SleepTime   = 10){
  
  QueueFinished <- F
  NrTrials      <- 0
  while((!QueueFinished) & NrTrials <= MaxNrTrials) {
    NrTrials <- NrTrials + 1
    cat("Checking queue, trial", NrTrials, "out of", MaxNrTrials, "\n")
    QueueStatus   <- system('squeue -u hzudohna',  intern = T)
    QueueFinished <- length(grep("batch", QueueStatus)) == 0
    Sys.sleep(SleepTime)
  }
  if (!QueueFinished){
    cat("Queue did not finish in", NrTrials*SleepTime, "seconds!\n")
    cat("Stopped checking\n")
  } else {
    cat("Queue finished!\n")
  }
  QueueFinished
}