# General description:
#
#    This function calculates the number of clipped bp from a cigar string

# Arguments:
#   
#    CigarString: sam file cigar string

# Output:
#   
#    number of bp clipped (from both sides)

# Comment:

NrClippedFromCigar <- function(CigarString){
  Nrs <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Nrs <- Nrs[c(1, length(Nrs))]
  Ltrs <- strsplit(CigarString, "[1-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  Ltrs <- Ltrs[c(1, length(Ltrs))]
  (Ltrs == "S") %*% Nrs
}



