# General description:
#
#    This function determines the aligned sequence of a read from the cigar 
#    string and the original sequence

# Arguments:
#   
#    CigarString: sam file cigar string
#    Seq: DNAstring set (created by scanBam)

# Output:
#   
#    SeqVector: aligned sequence as vector of upper case letters

# Comment:

SeqFromCigar <- function(CigarString, Seq){
  SeqV <- strsplit(as.character(Seq), "")[[1]]
  Nrs <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Ltrs <- strsplit(CigarString, "[0-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  StartV <- if(Ltrs[1] == "S") 2 else 1
  EndV   <- if(Ltrs[length(Ltrs)] == "S") length(Ltrs) - 1 else length(Ltrs)
  StartS <- if(Ltrs[1] == "S") Nrs[1] + 1 else 1
  SeqVector <- NULL
  for (i in StartV:EndV){
    Seq2Append <- switch(Ltrs[i], M = SeqV[StartS:(StartS + Nrs[i] - 1)],
                         D = rep("-", Nrs[i]),
                         I = NULL)
    SeqVector <- c(SeqVector, Seq2Append)
    StartS    <- StartS + Nrs[i] * (Ltrs[i] != "D")
  }
  SeqVector
}
