##############################################
#
# General description:
#
#   The following function obtains a numeric value of the length of 
#   a structural variant from a vcf INFO column
#   

# Input:
#
#     x: entry from INFO column in vcf file

# Output:
#   
#     numeric value for SV length

##############################################

GetFromVcfINFO_SVLength <- function(x){
  Split1 <- strsplit(x, ";")[[1]]
  LengthPart <- grep("SVLEN=", Split1, value = T)
  if (length(LengthPart) > 0){
    as.numeric(strsplit(LengthPart, "=")[[1]][2])
  } else {
    NA
  }
}