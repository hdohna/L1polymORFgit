# General description:
#
#    This function calculates the length of a mapped read from the cigar string

# Arguments:
#   
#    CigarString: sam file cigar string

# Output:
#   
#    number of bp clipped (from both sides)

# Comment:

ReadLengthFromCigar <- function(CigarString){
  Nrs <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Ltrs <- strsplit(CigarString, "[0-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  blnPlus  <- Ltrs %in% c("M", "D")
  sum(Nrs[blnPlus])
}