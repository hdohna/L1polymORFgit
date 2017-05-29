# General description:
#
#    This function determines the aligned phred score of a read from the cigar 
#    string and the original phred score sequence

# Arguments:
#   
#    CigarString: sam file cigar string
#    Phred: character string of phred score

# Output:
#   
#    PhredVector: aligned sequence as vector of upper case letters

# Comment:

PhredFromCigar <- function(CigarString, Phred){
  PhredV <- strsplit(as.character(Phred), "")[[1]]
  Nrs  <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Ltrs <- strsplit(CigarString, "[0-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  StartV <- if(Ltrs[1] == "S") 2 else 1
  EndV   <- if(Ltrs[length(Ltrs)] == "S") length(Ltrs) - 1 else length(Ltrs)
  StartS <- if(Ltrs[1] == "S") Nrs[1] + 1 else 1
  PhredVector <- NULL
  for (i in StartV:EndV){
    Phred2Append <- switch(Ltrs[i], M = PhredV[StartS:(StartS + Nrs[i] - 1)],
                         D = rep("-", Nrs[i]),
                         I = NULL)
    PhredVector <- c(PhredVector, Phred2Append)
    StartS    <- StartS + Nrs[i] * (Ltrs[i] != "D")
  }
  sapply(PhredVector, function(x) strtoi(charToRaw(x), base = 16L) - 33,
         USE.NAMES = F)
  
}
