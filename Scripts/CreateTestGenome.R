# The following script creates a reference genome to test LINE-1 detection
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)
#load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
load('/labs/dflev/hzudohna/RefSeqData/GRanges_L1_1000Genomes.RData')


# Read in L1HS consensus
#L1HSconsensus <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fas")
L1HSconsensus <- read.fasta("/labs/dflev/hzudohna/RefSeqData/Homo_sapiens_L1_consensus.fas")
L1HS_DNAString <- DNAString(paste(L1HSconsensus[[1]], collapse = ""))

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Get all chromosome number
AllChrs <- paste("chr", c(1:22, "X", "Y"), sep = "")


# Define function to create text for DNAString with insertions
CreateInsertTxt <- function(idx){
  Pos        <- L1_1000G$POS[idx]
  L1Starts   <- L1_1000G$L1StartNum[idx]
  L1Ends     <- L1_1000G$L1EndNum[idx]
  GenRanges  <- paste(c(1, Pos[-length(Pos)] + 1), Pos, sep = ":")
  L1Ranges   <- paste(L1Starts, L1Ends, sep = ":")
  L1Parts    <- paste(paste("L1HS_DNAString[", L1Ranges), "]")
  ChromParts <- paste(paste(CurrentChrom, "[", GenRanges), "]")
  NewDNASt_txt <- paste("c(", paste(paste(ChromParts, L1Parts, sep = ", "), collapse = ", "),
                        ")")
  
}

# Loop through the first 5 sample columns and generate reference genomes
# with the same insertions as the current
for (x in SampleColumns[6:26]){
  
  cat("******   Simulating genome", x, "   **********\n")
  
  # Path for genome fasta file path
  FastaPath <- "D:/L1polymORF/Data/hg19SampleGenome.fa"
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
  idxL1   <- which(L1_1000G[,x] > 0)
  ChrNrs <- L1_1000G$CHROM[idxL1]
  Chroms <- paste("chr", ChrNrs, sep = "")
  UniqueChroms <- unique(Chroms)

  # Loop over chromosomes and generate insertions
  for (Chr in AllChrs){
    cat("Processing", Chr, "\n")
    CurrentChrom  <- paste('BSgenome.Hsapiens.UCSC.hg19[["', Chr, '"]]', sep = "")
    NewDNASt_txt1 <- CurrentChrom
    NewDNASt_txt2 <- CurrentChrom
    if(Chr %in% UniqueChroms){
      idxChr     <- idxL1[Chr == Chroms]
      Count      <- L1_1000G[idxChr,x]
      bln2       <- Count == 2
      
      # Create indices for both homologous chromosome
      idx1       <- idxChr[bln2]
      idx2       <- idxChr[bln2]
      HaploSample <- sample(c(T, F), sum(!bln2), replace = T)
      idx1 <- c(idx1, idxChr[which(!bln2)[HaploSample]])
      idx2 <- c(idx2, idxChr[which(!bln2)[!HaploSample]])
      
      # Create insertion patterns for both homologous chromosome
      if (length(idx1) > 0) NewDNASt_txt1 <- CreateInsertTxt(idx1)
      if (length(idx2) > 0) NewDNASt_txt2 <- CreateInsertTxt(idx2)
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
  
}
