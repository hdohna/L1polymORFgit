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
  L1PosOrder <- order(L1_1000G$POS[idx])
  idx <- idx[L1PosOrder]
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

# Index of L1 in current genome x with L1 start and end values
idxL1   <- which(L1_1000G[,x] > 0 & 
                   (!is.na(L1_1000G$L1StartNum)) &
                   (!is.na(L1_1000G$L1EndNum)) )
ChrNrs <- L1_1000G$CHROM[idxL1]
Chroms <- paste("chr", ChrNrs, sep = "")
UniqueChroms <- unique(Chroms)
idxMismatchList <- list()
TotalL1Match    <- 0
TotalL1MisMatch <- 0
for (i in 1:length(idxList)){
  CurrentChrom <- names(idxList)[i]
  idxMismatchList[[i]] <- list()
  idxMismatchList[[i]][[1]] <- NULL
  idxMismatchList[[i]][[2]] <- NULL
  # Analyze haplotype 1
  cat("\n**********   Analyzing", CurrentChrom, "     ***********\n")
  ChrIdx <- which(AllChrs == CurrentChrom)
  if (length(idxList[[i]]$idx1) > 0) {
    Chr1   <- scan(FastaPath1, skip = 2*ChrIdx - 1, what = character(), nlines = 1)
    cat("Read haplotype 1 of length", nchar(Chr1), "\n")
    blnCheck1 <- CheckL1(idxList[[i]]$idx1, Chr1)
    cat("L1 checked in haplotype 1:", length(blnCheck1), "\n")
    cat("Correct L1 sequences:", sum(blnCheck1), "\n")
    cat("Incorrect L1 sequences:", sum(!blnCheck1), "\n")
    TotalL1Match    <- TotalL1Match + sum(blnCheck1)
    TotalL1MisMatch <- TotalL1MisMatch + sum(!blnCheck1)
    if(sum(!blnCheck1) > 0){
      idxMismatchList[[i]][[1]] <- idxList[[i]]$idx1[!blnCheck1]
      names(idxMismatchList[[i]][[1]]) <- "idx1"
    }
  } else {
    cat("No L1 in haplotype 1\n")
  }

  # Analyze haplotype 2
  if (length(idxList[[i]]$idx2) > 0) {
    Chr2 <- scan(FastaPath2, skip = 2*ChrIdx - 1, what = character(), nlines = 1)
    cat("Read haplotype 2 of length", nchar(Chr2), "\n")
    blnCheck2 <- CheckL1(idxList[[i]]$idx2, Chr2)
    cat("L1 checked in haplotype 2:", length(blnCheck2), "\n")
    cat("Correct L1 sequences:", sum(blnCheck2), "\n")
    cat("Incorrect L1 sequences:", sum(!blnCheck2), "\n")
    TotalL1Match    <- TotalL1Match + sum(blnCheck2)
    TotalL1MisMatch <- TotalL1MisMatch + sum(!blnCheck2)
    if(sum(!blnCheck2) > 0){
      idxMismatchList[[i]][[2]] <- idxList[[i]]$idx2[!blnCheck2]
      names(idxMismatchList[[i]][[2]]) <- "idx2"
    }
  } else {
    cat("No L1 in haplotype 2\n")
  }
}
# names(idxMismatchList) <- names(idxList)
# rm(list = c("Chr1", "Chr2"))
cat("\n**************************************************************\n")
cat("*                                                            *\n")
cat("*   Total match:", TotalL1Match, "                           *\n")
cat("*   Total mismatch:", TotalL1MisMatch, "                       *\n")
cat("*                                                            *\n")
cat("**************************************************************\n")



# Save idxList
# OutPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/L1MismatchInfo_", 
#                  x, ".RData", sep = "")
# cat("Saving imismatch to", OutPath, " ... ")
# save(list = "idxList", file = OutPath)
# cat("done!")

