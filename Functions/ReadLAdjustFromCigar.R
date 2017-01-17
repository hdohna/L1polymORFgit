# General description:
#
#    This function calculates the number of bp from a cigar string that should
#    be used to turn sequence length into length of mapped read

# Arguments:
#   
#    CigarString: sam file cigar string

# Output:
#   
#    number of bp clipped (from both sides)

# Comment:

ReadLAdjustFromCigar <- function(CigarString){
  Nrs <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Ltrs <- strsplit(CigarString, "[1-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  ULtrs <- unique(Ltrs)
  LtrSigns <- rep(0, length(ULtrs))
  LtrSigns[LtrSigns %in% c("S", "I")] <- -1
  sapply(unique(Ltrs), function(x) (Ltrs))
}



