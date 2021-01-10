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


SummarizeQueue <- function(){
  
    QueueStatus <- system('squeue -u hb54',  intern = T)
    idxBatch    <- grep("batch", QueueStatus)
    idxInteract <- grep("interact", QueueStatus)
    idxRunning  <- grep(" R ", QueueStatus)
    idxPending  <- grep(" PD ", QueueStatus)
    cat("\n\n******   Batch jobs:   **********\n")
    cat(length(intersect(idxBatch, idxRunning)), "jobs running and ")
    cat(length(intersect(idxBatch, idxPending)), "jobs pending\n\n")
    cat("******   Interactive jobs:   **********\n")
    cat(length(intersect(idxInteract, idxRunning)), "jobs running and ")
    cat(length(intersect(idxInteract, idxPending)), "jobs pending\n")
}