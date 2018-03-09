##############################################
#
# General description:
#
#   The following function gets the columns from a sngleton
#   file

# Input:
#
#     SingletonFile: File as created by Singleton_awk_cmds

# Output:
#   


# Comments:
#    

##############################################


GetSingletonColumns <- function(SingletonFile){
  
  # Subtract first row
  minusFirstRow <- SingletonFile$V10 - 20031
  if (SingletonFile$V1[1] == 1){
    beyondRow  <- minusFirstRow %/% 10016 >= 5226911
    PredictCol <- ceiling(((minusFirstRow - beyondRow*11401) %% 10016)/4)
  } else {
    PredictCol <- ceiling((minusFirstRow %% 10016)/4)
  }
  PredictCol
}