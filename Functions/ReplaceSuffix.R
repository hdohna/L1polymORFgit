##############################################
#
# General description:
#
#   The following function replaces the suffix of a file name (indicated by .)

# Input:
#
#     FileName: path to file
#     Replacement: new suffix 

# Output (as list):
#   
#     path to a new file with replaced suffix

# Comments:

##############################################

ReplaceSuffix <- function(FileName, Replacement){
  NameSplit <- strsplit(FileName, "\\.")[[1]]
  paste(c(NameSplit[-length(NameSplit)], Replacement), collapse = ".")
}
