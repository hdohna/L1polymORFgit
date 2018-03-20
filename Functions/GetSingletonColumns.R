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
  minusRow <- SingletonFile$V10 - 20031
  
  # Perform correction for chromosome 1 (it contains 11401 characters in row 5226911)
  if (SingletonFile$V1[1] == 1){
    beyondRow  <- minusRow %/% 10016 >= 5226911
    nrChar     <- (minusRow - beyondRow * 11401) %% 10016
  } else {
    nrChar     <- minusRow %% 10016
  }
  Allele     <- nrChar %% 4
  PredictCol <- ceiling(nrChar/4)
  data.frame(Col = PredictCol, Allele = Allele)
}