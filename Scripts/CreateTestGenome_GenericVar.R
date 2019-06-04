# The following script is a generic script to create a reference 
# genome to test LINE-1 detection

# GenericID has to be replaced to run the script
x <- "GenericID"

library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

load('/labs/dflev/hzudohna/RefSeqData/GRanges_L1_1000Genomes.RData')
load('/labs/dflev/hzudohna/RefSeqData/ChromLengthsHg19.Rdata')

# Read in L1HS consensus
#L1HSconsensus <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas")
L1HSconsensus <- read.fasta("/labs/dflev/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas")
L1HS_DNAString <- DNAString(paste(L1HSconsensus[[1]], collapse = ""))
L1Length <- length(L1HS_DNAString)

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Create a column for proportion different from consensus
L1_1000G$PropDiff   <- runif(nrow(L1_1000G), max = 0.2)

# Get all chromosome number
AllChrs <- paste("chr", c(1:22, "X", "Y"), sep = "")


# Define function to create text for DNAString with insertions
CreateInsertTxt <- function(idx, ChrLength){
  Pos        <- L1_1000G$POS[idx]
  L1Starts   <- L1_1000G$L1StartNum[idx]
  L1Ends     <- L1_1000G$L1EndNum[idx]
  PosOrder   <- order(Pos)
  Pos        <- Pos[PosOrder]
  L1Starts   <- L1Starts[PosOrder]
  L1Ends     <- L1Ends[PosOrder]
  GenRanges  <- paste(c(1, Pos + 1), c(Pos, ChrLength), sep = ":")
  L1Ranges   <- paste(L1Starts, L1Ends, sep = ":")
  L1Parts    <- paste(paste("L1HS_DNAString[", L1Ranges), "]")
  ChromParts <- paste(paste(CurrentChrom, "[", GenRanges), "]")
  NewDNASt_txt <- paste("c(", paste(paste(ChromParts, L1Parts, sep = ", "), 
                                    collapse = ", "),
                        ")")
}


# Create a list of L1HS with specified deviations from consensus
cat("Creating variable L1HS ...")
L1HSList <- lapply(L1_1000G$PropDiff, function(x){
    NrDiff  <- round(x * L1Length);
    DiffPos <- sample.int(L1Length, NrDiff);
    NewNuc  <- sample(c("A", "C", "G", "T"), size = NrDiff,
                      replace = T);
    NewL1Hs <- L1HS_DNAString;
    NewL1Hs[DiffPos] <- DNAString(paste(NewNuc, collapse = ""));
    STCompare     <- compareStrings(NewL1Hs[DiffPos], L1HS_DNAString[DiffPos]);
    STCompareVect <- s2c(STCompare);
    idxSame       <- which(STCompareVect != "?")
    for (Nuc in c("A", "C", "G", "T")){
      blnRepl <- STCompareVect == Nuc
      ReplaceNucs <- sample(setdiff(c("A", "C", "G", "T"), Nuc), sum(blnRepl),
                            replace = T)
      NewL1Hs[DiffPos][blnRepl] <- DNAString(paste(ReplaceNucs, collapse = ""))

    }
  NewL1Hs
})
cat("done!\n")

# Define function to create text for DNAString with insertions
# with deviations from consensus
CreateInsertTxtVar <- function(idx){
  Pos        <- L1_1000G$POS[idx]
  GenRanges  <- paste(c(1, Pos[-length(Pos)] + 1), Pos, sep = ":")
  L1Parts    <- paste(paste("L1HSList[[", idx), "]]")
  ChromParts <- paste(paste(CurrentChrom, "[", GenRanges), "]")
  NewDNASt_txt <- paste("c(", paste(paste(ChromParts, L1Parts, sep = ", "), collapse = ", "),
                        ")")
}

# Loop through the first sample columns and generate reference genomes
# with the same insertions as the current
idxList <- list()

  cat("******   Simulating genome", x, "   **********\n")
  ListNames <- NULL

  # Path for genome fasta file path
  FastaPath1 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19_",
                     x, "_Haplo1.fa", sep = "")
  FastaPath2 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19_",
                      x, "_Haplo2.fa", sep = "")
  if (file.exists(FastaPath1)){
    file.remove(FastaPath1)
  }
  if (file.exists(FastaPath2)){
    file.remove(FastaPath2)
  }
  
  # Open connections for fasta file to write chromosomes
  GenCon1 <- file(FastaPath1, open = "a")
  GenCon2 <- file(FastaPath2, open = "a")
  
  # Index of L1 in current genome and all chromosomes with L1
  idxL1   <- which(L1_1000G[,x] > 0 & 
                     (!is.na(L1_1000G$L1StartNum)) &
                     (!is.na(L1_1000G$L1EndNum)) )
  ChrNrs <- L1_1000G$CHROM[idxL1]
  Chroms <- paste("chr", ChrNrs, sep = "")
  UniqueChroms <- unique(Chroms)

  # Loop over chromosomes and generate insertions
  idxLevel2 <- 0
  for (Chr in AllChrs){
    cat("Processing", Chr, "\n")
    idxLevel2 <- idxLevel2 + 1
    CurrentChrom  <- paste('BSgenome.Hsapiens.UCSC.hg19[["', Chr, '"]]', sep = "")
    NewDNASt_txt1 <- CurrentChrom
    NewDNASt_txt2 <- CurrentChrom
    if(Chr %in% UniqueChroms){
      idxChr     <- idxL1[Chr == Chroms]
      Count      <- L1_1000G[idxChr,x]
      bln2       <- Count == 2
      
      # Create indices for both homologous chromosomes
      idx1       <- idxChr[bln2]
      idx2       <- idxChr[bln2]
      HaploSample <- sample(c(T, F), sum(!bln2), replace = T)
      idx1 <- c(idx1, idxChr[which(!bln2)[HaploSample]])
      idx2 <- c(idx2, idxChr[which(!bln2)[!HaploSample]])
      
      # Add indices of haplotypes to list that keeps track of them
      idxList[[idxLevel2]] <- list(idx1 = idx1, idx2 = idx2)
      ListNames <- c(ListNames, Chr)
      
      # # Create insertion patterns for both homologous chromosomes
      # if (length(idx1) > 0) NewDNASt_txt1 <- CreateInsertTxt(idx1, 
      #                                                        ChromLengthsHg19[Chr])
      # if (length(idx2) > 0) NewDNASt_txt2 <- CreateInsertTxt(idx2, 
      #                                                        ChromLengthsHg19[Chr])
      # Create insertion patterns for both homologous chromosome
      if (length(idx1) > 0) NewDNASt_txt1 <- CreateInsertTxtVar(idx1)
      if (length(idx2) > 0) NewDNASt_txt2 <- CreateInsertTxtVar(idx2)
    }
    NewDNASt1     <- eval(parse(text = NewDNASt_txt1))
    NewDNASt2     <- eval(parse(text = NewDNASt_txt2))
    NewDNASt_char1 <- as.character(NewDNASt1)
    NewDNASt_char2 <- as.character(NewDNASt2)
    ChromName <- paste(">", substr(Chr, 4, nchar(Chr)), sep = "")
    writeLines(text = c(ChromName, NewDNASt_char1), con = GenCon1)
    writeLines(text = c(ChromName, NewDNASt_char2), con = GenCon2)
  }
  close(GenCon1)  
  close(GenCon2)  
  names(idxList) <- ListNames
  
# Save idxList
OutPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/HaploIdxListVar_", 
                 x, ".RData", sep = "")
save(list = "idxList", file = OutPath)

