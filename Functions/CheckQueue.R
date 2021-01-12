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
                       SleepTime   = 10,
                       JobIDs = NULL){
  
  QueueFinished <- F
  NrTrials      <- 0
  while((!QueueFinished) & NrTrials < MaxNrTrials) {
    
    NrTrials <- NrTrials + 1
    cat("Checking queue, trial", NrTrials, "out of", MaxNrTrials, 
        ", ", round((NrTrials - 1) * SleepTime / 60, 2), "min elapsed - ")
    QueueStatus   <- system('squeue -u hb54',  intern = T)
    if(!is.null(JobIDs)){
      QueueStatus <- sapply(JobIDs, function(x) grep(x, QueueStatus, value = T))
    }
    idxRunning    <- grep(" R ", QueueStatus)
    idxPending    <- grep(" PD ", QueueStatus)
    QueueFinished <- (length(grep("bash", QueueStatus)) + 1) == length(QueueStatus)
    if (!QueueFinished) cat("queue not yet finished!", 
                            length(idxRunning),
                            "jobs are running and",
                            length(idxPending),
                            "jobs are pending.\n")
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