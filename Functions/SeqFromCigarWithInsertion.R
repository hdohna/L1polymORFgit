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
#    InsMat: matrix of insertion start and end

# Comment:
#
#    Similar to function SeqFromCigar but keeps track of insertion positions

SeqFromCigarWithInsertion <- function(CigarString, Seq){

  SeqV <- strsplit(as.character(Seq), "")[[1]]
  Nrs  <- as.numeric(strsplit(CigarString, "[A-Z]")[[1]])
  Ltrs <- strsplit(CigarString, "[0-9]")[[1]]
  Ltrs <- Ltrs[Ltrs != ""]
  
  # Create objects for skipping Ds or Is
  NrsWithoutD       <- Nrs
  NrsWithoutI       <- Nrs
  blnD              <- Ltrs == "D"
  blnI              <- Ltrs == "I"
  NrsWithoutD[blnD] <- 0
  NrsWithoutI[blnI] <- 0
  
  # Create vectors of start and end indices on sequence vector
  idxSeqStart <- cumsum(c(0, NrsWithoutD[-length(NrsWithoutD)])) + 1
  idxSeqEnd   <- cumsum(NrsWithoutD)

  # Create vectors of start and end indices on genome
  idxConsStart <- cumsum(c(0, NrsWithoutI[-length(NrsWithoutI)])) + 1
  idxConsEnd   <- cumsum(NrsWithoutI)
  
  # Get indices along the sequences
  if(Ltrs[1] == "S") {
    idxSeqStart  <- idxSeqStart[-1]
    idxSeqEnd    <- idxSeqEnd[-1]
    
    # Indices along the consensus sequence
    idxConsStart <- idxConsStart[-1]
    idxConsEnd   <- idxConsEnd[-1]
    Nrs          <- Nrs[-1]
    blnD         <- blnD[-1]
    blnI         <- blnI[-1]
  }
  if(Ltrs[length(Ltrs)] == "S"){
    Li <- length(idxSeqStart)
    idxSeqStart  <- idxSeqStart[-Li]
    idxSeqEnd    <- idxSeqEnd[-Li]
    
    # Indices along the consensus sequence
    idxConsStart <- idxConsStart[-Li]
    idxConsEnd   <- idxConsEnd[-Li]
    Nrs          <- Nrs[-Li]
    blnD         <- blnD[-Li]
    blnI         <- blnI[-Li]
    
  }
  ConsDiff     <- idxConsStart[1] - 1
  idxConsStart <- idxConsStart - ConsDiff
  idxConsEnd   <- idxConsEnd   - ConsDiff
  
  # Indices  to build sequence
  idxM   <- which(!(blnD | blnI))
  
  # Loop to build sequences
  SeqVector <- NULL
  for (i in idxM){
    SeqVector <- c(SeqVector, SeqV[idxSeqStart[i]:idxSeqEnd[i]])
  }

  # Loop and add blanks
  for (i in which(blnD)){
    SeqVector <- append(SeqVector, rep("-", Nrs[i]), idxConsStart[i] - 1)
  }
  SeqVector
}

