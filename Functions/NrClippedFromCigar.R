# General description:
#
#    This function calculates the number of clipped bp from a cigar string

# Arguments:
#   
#    CigarString: sam file cigar string

# Output:
#   
#    vector (length 2) of number of left and right bp clipped 

# Comment:

NrClippedFromCigar <- function(CigarString){
  Nrs <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Nrs <- Nrs[c(1, length(Nrs))]
  Ltrs <- strsplit(CigarString, "[0-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  Ltrs <- Ltrs[c(1, length(Ltrs))]
  (Ltrs == "S") * Nrs
}



