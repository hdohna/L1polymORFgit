##############################################
#
# General description:
#
#   The following function gets L1 sequences, genomic coordinates and flanking
#   sequences from a vector of accession numbers (or a list of sequences)

# Input:
#
#  Folders and file names:
#     L1Consens:     L1HS consensus sequence
#     AccNrs:        vector with accession numbers to get sequences
#     Seqs:          DNAStringSet of clone sequences (if it is null, sequences will be 
#                    obtained using AccNrs)
#     MinMatchWidth: integer giving the minimum width of a match to L1 sequence
#     FlankSize:     integer giving the length of sequences flanking L1 insertions 
#                    [bp] used to locate insertions


# Output:
#   
#    OutDF: data frame containing the following colums
#        L1Seq         : L1 sequence
#        L1SeqFlank5p  : 5' flanking sequence
#        L1SeqFlank3p  : 3' flanking sequence
#        Strand        : strand


##############################################

getNonRefL1 <- function(L1Consens,
  AccNrs, Seqs = NULL, MinMatchWidth = 5500, FlankSize = 100,
  Chromosomes, RefGenome = BSgenome.Hsapiens.UCSC.hg38,
  blnLocateL1inRef = F){
  
  cat("\n\n*******************************************************\n")
  cat("**                                                   **\n")
  cat("**    Running function getNonRefL1 ...               **\n")
  cat("**                                                   **\n")
  cat("*******************************************************\n\n")
  
  # Check that chromosomes and accession number have the the same length
  if (length(AccNrs) != length(Chromosomes)){
    stop("AccNrs and Chromosomes have to match!")
  }
  
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
  cat("\n***** Searching for match of consensus L1 in clone *****\n\n")
  SeqsRV <- reverseComplement(Seqs)
  
  # Initialize output objects
  OutDF <- data.frame(
    L1Seq             = rep(NA, length(Seqs)),
    L1SeqFlank5p      = rep(NA, length(Seqs)),
    L1SeqFlank3p      = rep(NA, length(Seqs)),
    L1SeqFlank5p2x    = rep(NA, length(Seqs)),
    L1SeqFlank3p2x    = rep(NA, length(Seqs)),
    Strand            = rep(NA, length(Seqs)),
    start_Clone       = rep(NA, length(Seqs)),
    end_Clone         = rep(NA, length(Seqs)),
    start_Ref         = rep(NA, length(Seqs)),
    end_Ref           = rep(NA, length(Seqs)),
    FlankStart        = rep(NA, length(Seqs)),
    InsertionStartRel = rep(NA, length(Seqs)),
    NrNuc5p           = rep(NA, length(Seqs)),
    NrNuc3p           = rep(NA, length(Seqs))
  )
  
  # Loop through sequences that have no coordinates in the reference genome
  # and get the sequences through pairwise alignment
  for (i in 1:length(Seqs)){
    cat("Processing", names(Seqs)[i], "sequence", i, "out of", length(Seqs), "\n")
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
      Seq5P2x <- SubjectSeq[[1]][max(1, (S - 2*FlankSize)):max(1, S - FlankSize - 1)] 
      Seq3P2x <- SubjectSeq[[1]][min(length(SubjectSeq[[1]]), E + 1 + FlankSize):
                                        min(length(SubjectSeq[[1]]), (E + 2*FlankSize))] 
      OutDF$L1SeqFlank5p[i]    <- as.character(Seq5P)
      OutDF$L1SeqFlank3p[i]    <- as.character(Seq3P)
      OutDF$L1SeqFlank5p2x[i]  <- as.character(Seq5P2x)
      OutDF$L1SeqFlank3p2x[i]  <- as.character(Seq3P2x)
      OutDF$start_Clone[i]     <- S
      OutDF$end_Clone[i]       <- E
    } else {
      cat("L1 consensus could not be matched to clone", 
          AccNrs[i], "\n")
    }
  }
  
  # Looking for insertion site in reference
  if (blnLocateL1inRef){
    cat("\n***** Looking for insertion site in reference ... *****\n\n")
    
    # Get the total length of the flanking sequence
    FlankTotal <- nchar(OutDF$L1SeqFlank5p2x) +
      nchar(OutDF$L1SeqFlank3p2x)
    
    # Loop through flanking sequences above minimal length, locate them on the
    # reference genome and get insert location descriptors
    for (i in which(FlankTotal > FlankSize)){
      cat("Analyzing insertion", i, "of", nrow(OutDF), "\n")
      if (OutDF$Strand[i] == "+") {
        LPattern <- OutDF$L1SeqFlank5p2x[i]
        RPattern <- OutDF$L1SeqFlank3p2x[i]
        InsertionSite <- paste(OutDF$L1SeqFlank5p[i],
                               OutDF$L1SeqFlank3p[i], sep = "")
      } else {
        DNASt5P   <- DNAString(OutDF$L1SeqFlank5p[i])
        DNASt3P   <- DNAString(OutDF$L1SeqFlank3p[i])
        DNASt5P2x <- DNAString(OutDF$L1SeqFlank5p2x[i])
        DNASt3P2x <- DNAString(OutDF$L1SeqFlank3p2x[i])
        LPattern  <- reverseComplement(DNASt3P2x)
        RPattern  <- reverseComplement(DNASt5P2x)
        InsertionSite <- paste(reverseComplement(DNASt3P),
                               reverseComplement(DNASt5P), sep = "")
      }
      Chrom    <- Chromosomes[i]
      ChromSeq <- RefGenome[[Chrom]]
      String   <- matchLRPatterns(LPattern, RPattern, 
                                  max.gaplength = 4 * FlankSize, ChromSeq,
                                  max.Lmismatch = 5, max.Rmismatch = 5)
      InsertStart <- NA
      if (length(String) > 0){
        c("Matched L1 flanks from clone to reference\n")
        OutDF$FlankStart[i]        <- start(String)
        OutDF$InsertionStartRel[i] <- 2 * FlankSize
        OutDF$start_Ref[i] <- start(String) + 2 * FlankSize
        OutDF$end_Ref[i]   <- start(String) + 2 * FlankSize + 1
        OutDF$NrNuc5p[i]           <- 0
        OutDF$NrNuc3p[i]           <- 0
        pwA <- pairwiseAlignment(InsertionSite, String, type = "local")
        Indel <- subsetByOverlaps(unlist(pwA@subject@indel), 
                                  IRanges(start = FlankSize,
                                          end = FlankSize))
        if (length(Indel) > 0){
          OutDF$InsertionStartRel[i] <- FlankSize + start(Indel) - 1
          OutDF$start_Ref[i]         <- start(String) + FlankSize + start(Indel) - 1
          OutDF$end_Ref[i]         <- start(String) + FlankSize + start(Indel)
          OutDF$NrNuc5p[i]           <- FlankSize - start(Indel)
          OutDF$NrNuc3p[i]           <- end(Indel) - FlankSize
          
        }
      }
      
    } # End of for (i in which(FlankTotal > FlankSize))
  } # End of if (blnLocateL1inRef)
  
  # Return output data.frame
  return(OutDF)
  
}

