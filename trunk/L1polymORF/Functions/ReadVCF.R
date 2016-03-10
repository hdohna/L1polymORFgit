##############################################
#
# General description:
#
#   The following function reads a vcf file as table
#   

# Input:
#
#     VcfFile: character string providing path to vcf file

# Output:
#   
#     table with content of vcf file

##############################################

# Function to read vcf file as delimited file
ReadVCF <- function(VcfFile) {
  VCFLines <- readLines(VcfFile)
  read.delim(file = VcfFile, skip = max(grep("##", VCFLines)), as.is = T)
}
