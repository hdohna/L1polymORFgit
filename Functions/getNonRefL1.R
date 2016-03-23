##############################################
#
# General description:
#
#   The following function gets L1 sequences, genomic coordinates and flanking
#   sequences from a vector of accession numbers (or a list of sequences)

# Input:
#
#  Folders and file names:
#     L1Consens: L1HS consensus sequence
#     AccNrs:    vector with accession numbers to get sequences
#     Seqs:      DNAStringSet of clone sequences (if it is null, sequences will be 
#                obtained using AccNrs)
#     MinMatchWidth: integer giving the minimum width of a match to L1 sequence


# Output:
#   
#    OutDF: data frame containing the following colums
#        L1Seq         : L1 sequence
#        L1SeqFlank5p  : 5' flanking sequence
#        L1SeqFlank3p  : 3' flanking sequence
#        Strand        : strand


##############################################

getNonRefL1 <- function(L1Consens,
  AccNrs, Seqs = NULL, MinMatchWidth = 5500){
  
  # Get sequences from Genbank if no sequences are provided
  if (is.null(Seqs)){
    cat("***** Downloading sequences from Genbank ... \n")
    choosebank("genbank")
    Seqs <- sapply(AccNrs, function(AccNr){
      cat("Processing", AccNr, "\n")
      x <- query(listname = "L1", paste("AC=", AccNr, sep = ""))
      Seq <- getSequence(x$req)
      paste(Seq[[1]], collapse = "")
    })
    closebank()
    names(Seqs) <- AccNrs
    cat("Done!  *****\n")
    Seqs <- toupper(Seqs)
    
    # Find sequences that are matched by flanks
    Seqs   <- DNAStringSet(Seqs)
  }

  # Check consistency of input
  if (length(Seqs) != length(AccNrs)){
    stop("Input vectors Seqs and AccNrs have to match!\n")
  }
  
  # Match each L1HS to a locus 
  cat("***** Searching for match of consensus L1  *****\n")
  SeqsRV <- reverseComplement(Seqs)
  
  # Initialize output objects
  OutDF <- data.frame(
    L1Seq         = rep(NA, length(Seqs)),
    L1SeqFlank5p  = rep(NA, length(Seqs)),
    L1SeqFlank3p  = rep(NA, length(Seqs)),
    Strand        = rep(NA, length(Seqs))
  )

  
  # Loop through sequences that have no coordinates in the reference genome
  # and get the sequences through pairwise alignment
  for (i in 1:length(Seqs)){
    cat("Processing sequence", i, "out of", length(Seqs), "\n")
    SubjectSeq <- Seqs[i]
    lpA <- pairwiseAlignment(L1Consens, SubjectSeq, type = "local")
    if (width(lpA@subject) > MinMatchWidth){
      Strand <- "+"
    } else {
      SubjectSeq <- SeqsRV[i]
      lpA        <- pairwiseAlignment(L1Consens, SubjectSeq, type = "local")
      Strand     <- "-"
    }
    if (width(lpA@subject) > MinMatchWidth){
      OutDF$L1Seq[i]    <- as.character(lpA@subject)
      OutDF$Strand[i]   <- Strand
      S <- start(lpA@subject@range)
      E <- end(lpA@subject@range)
      Seq5P <- SubjectSeq[[1]][max(1, (S - FlankSize)):(S - 1)] 
      Seq3P <- SubjectSeq[[1]][(E + 1):min(length(SubjectSeq[[1]]), (E + FlankSize))] 
      OutDF$L1SeqFlank5p[i]  <- as.character(Seq5P)
      OutDF$L1SeqFlank3p[i]  <- as.character(Seq3P)
    } else {
      cat("L1 consensus could not be matched to clone", 
          AccNrs[i], "\n")
    }
  }
  
  # Return output data.frame
  return(OutDF)
  
}

