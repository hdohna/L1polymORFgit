##############################################
#
# General description:
#
#   The following function extracts strand of an LINE1 insertion from 
#   information provided in 1000 Genome 
#   

# Input:
#
#     InfoString: character string providing info about LINE1


# Output:
#   
#     Strand ("-" or "+")

##############################################

GetStrandFrom1000GVcf <- function(InfoString){
  Split1 <- strsplit(InfoString, ";")[[1]]
  y <- grep("MEINFO=", Split1, value = T)
  if(length(y) > 0){
    substr(y, nchar(y), nchar(y))
  } else {
    "*"
  }
  
}