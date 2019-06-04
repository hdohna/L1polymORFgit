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
L1HSconsCat   <- paste(toupper(L1HSconsensus[[1]]), collapse = "")
L1HS_DNAString <- DNAString(L1HSconsCat)
L1Length <- length(L1HS_DNAString)

# Add numeric columns for L1 start and end
L1_1000G$L1StartNum <- as.numeric(as.character(L1_1000G$L1Start))
L1_1000G$L1EndNum   <- as.numeric(as.character(L1_1000G$L1End))

# Get all chromosome number
AllChrs <- paste("chr", c(1:22, "X", "Y"), sep = "")

# Function to Check whether L1 along the chromosome
CheckL1 <- function(idx, ChrSeq){
  L1Start   <- L1_1000G$L1StartNum[idx]
  L1End     <- L1_1000G$L1EndNum[idx]
  L1Widths  <- c(0, L1End - L1Start + 1)
  L1CumSum  <- cumsum(L1Widths)
  ChrStarts <- L1_1000G$POS[idx] + L1CumSum[-length(L1CumSum)] + 1
  ChrEnds   <- L1_1000G$POS[idx] + L1CumSum[-1]
  sapply(1:length(ChrStarts), function(i){
    substr(ChrSeq, ChrStarts[i], ChrEnds[i]) ==
    substr(L1HSconsCat, L1Start[i], L1End[i])
    
  })
}


# Path for genome fasta file path
FastaPath1 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19_",
                    x, "_Haplo1.fa", sep = "")
FastaPath2 <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/hg19_",
                    x, "_Haplo2.fa", sep = "")

# Path to list with indices
ListPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/HaploIdxList_", 
                 x, ".RData", sep = "")
load(ListPath)

cat("******   Checking genome", x, "   **********\n")

# Index of L1 in current genome x all chromosomes with L1
idxL1   <- which(L1_1000G[,x] > 0 & 
                   (!is.na(L1_1000G$L1StartNum)) &
                   (!is.na(L1_1000G$L1EndNum)) )
ChrNrs <- L1_1000G$CHROM[idxL1]
Chroms <- paste("chr", ChrNrs, sep = "")
UniqueChroms <- unique(Chroms)
length(idxList)
for (i in 1:length(idxList)){
  CurrentChrom <- names(idxList)[i]
  cat("\n**********   Analyzing", CurrentChrom, "     ***********\n")
  ChrIdx <- which(AllChrs == CurrentChrom)
  Chr1   <- scan(FastaPath1, skip = 2*ChrIdx - 1, what = character(), nlines = 1)
  blnCheck1 <- CheckL1(idxList[[i]]$idx1, Chr1)
  cat("Checked", length(blnCheck1), "L1 in haplotype 1\n")
  if (all(blnCheck1)){
    cat("All L1 have correct sequences!\n")
  } else {
    cat("Warning:", sum(!blnCheck1), "L1 sequences not correct!\n")
  }
  Chr2 <- scan(FastaPath2, skip = 2*ChrIdx - 1, what = character(), nlines = 1)
  blnCheck2 <- CheckL1(idxList[[i]]$idx2, Chr2)
  cat("Checked", length(blnCheck2), "L1 in haplotype 2\n")
  if (all(blnCheck2)){
    cat("All L1 have correct sequences!\n")
  } else {
    cat("Warning:", sum(!blnCheck2), "L1 sequences not correct!\n")
  }
}
cat("done!")

