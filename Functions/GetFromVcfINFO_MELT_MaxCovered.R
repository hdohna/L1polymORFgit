##############################################
#
# General description:
#
#   The following function obtains a numeric value of the last position covered 
#   by reads from a vcf INFO column
#   generated by MELT
#   

# Input:
#
#     x: entry from INFO column in vcf file generated by MELT

# Output:
#   
#     numeric value of proportion of a mobile element covered by reads

##############################################

GetFromVcfINFO_MELT_MaxCovered <- function(x){
  Split1 <- strsplit(x, ";")[[1]]
  DiffPart1 <- grep("DIFF=", Split1, value = T)
  if (length(DiffPart1) > 0){
    DiffPart2 <- strsplit(DiffPart1, ":")[[1]]
    DiffPart3 <- strsplit(DiffPart2[2], ",")[[1]]
    DiffPart4 <- strsplit(DiffPart3[length(DiffPart3)], "-")[[1]][1]
    as.numeric(substr(DiffPart4, 2, nchar(DiffPart4)))
  } else {
    NA
  }
}
