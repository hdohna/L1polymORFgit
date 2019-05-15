##############################################
#
# General description:
#
#   The following function obtains a numeric value of the genotype 
#   from a vcf INFO column
#   

# Input:
#
#     x: entry from genotype column in vcf file
#     GenoSplit: sign for splitting alleles ("/" for unphased and "|" for 
#                phased)

# Output:
#   
#     numeric value for genotype (0, 1, or 2)

##############################################

GetFromVcfGeno_GenoNum <- function(x, GenoSplit = "/"){
  Split1 <- strsplit(x, GenoSplit)[[1]]
  Split2 <- strsplit(Split1[2], ":")[[1]][1]
  sum(as.numeric(c(Split1[1], Split2)))
}