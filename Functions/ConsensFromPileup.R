# General description:
#
#    This function determines the consensus sequence from a pileup

# Arguments:
#   
#    Pileup: object created by the function pileup {Rsamtools}
#    NullSeq: character vector (dashes)

# Output:
#   
#    NullSeq: Character vector with dashes replaced by consensus

# Comment:

ConsensFromPileup <- function(Pileup, NullSeq = rep("-", max(Pileup$pos))){
  
  if (length(Pileup$pos) > 0){
    NucPosTab <- table(Pileup$pos, Pileup$nucleotide)
    NucPosTabNonZero <- NucPosTab > 0
    NrDiffNucs <- rowSums(NucPosTabNonZero)
    if (any(NrDiffNucs > 1)){
      cat("Ambiguous nucleotides at positions", 
          rownames(NucPosTab)[NrDiffNucs > 1], "\n")
    }
    ReplaceNuc <- tolower(colnames(NucPosTab))
    Consens <- apply(NucPosTab, 1, FUN = function(x) ReplaceNuc[which.max(x)])
    NullSeq[as.numeric(names(Consens))] <- Consens
    
  }
  NullSeq
}



