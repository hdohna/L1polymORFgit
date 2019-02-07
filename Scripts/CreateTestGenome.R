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


# Loop through the first 5 sample columns and generate reference genomes
# with the same insertions as the current
for (x in SampleColumns[1:5]){
  
  cat("******   Simulating genome", x, "   **********\n")
  
  # Path for genome fasta file path
  FastaPath <- "D:/L1polymORF/Data/hg19SampleGenome.fa"
  FastaPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19_",
                     x, ".fa", sep = "")
  if (file.exists(FastaPath)){
    file.remove(FastaPath)
  }
  GenCon <- file(FastaPath, open = "a")
  
  idxL1 <- which(L1_1000G[,x] > 0)
  Genome1 <- BSgenome.Hsapiens.UCSC.hg19
  ChrNrs <- L1_1000G$CHROM[idxL1]
  Chroms <- paste("chr", ChrNrs, sep = "")
  UniqueChroms <- unique(Chroms)
  Chr <- UniqueChroms[1]
  for (Chr in AllChrs){
    cat("Processing", Chr, "\n")
    CurrentChrom <- paste('BSgenome.Hsapiens.UCSC.hg19[["', Chr, '"]]', sep = "")
    NewDNASt_txt <- CurrentChrom
    if(Chr %in% UniqueChroms){
      idxChr    <- idxL1[Chr == Chroms]
      Pos       <- L1_1000G$POS[idxChr]
      L1Starts  <- L1_1000G$L1StartNum[idxChr]
      L1Ends    <- L1_1000G$L1EndNum[idxChr]
      GenRanges <- paste(c(1, Pos[-length(Pos)] + 1), Pos, sep = ":")
      L1Ranges <- paste(L1Starts, L1Ends, sep = ":")
      L1Parts  <- paste(paste("L1HS_DNAString[", L1Ranges), "]")
      ChromParts <- paste(paste(CurrentChrom, "[", GenRanges), "]")
      NewDNASt_txt <- paste("c(", paste(paste(ChromParts, L1Parts, sep = ", "), collapse = ", "),
                            ")")
    }
    NewDNASt <- eval(parse(text = NewDNASt_txt))
    NewDNASt_char <- as.character(NewDNASt)
    ChromName <- paste(">", substr(Chr, 4, nchar(Chr)), sep = "")
    writeLines(text = c(ChromName, NewDNASt_char), con = GenCon)
  }
  close(GenCon)  
  
}

