###############################################
#                                             #
#       Create alignments of reads            #
#         covering full-length L1             #
#                                             #
###############################################

# Get indices of insertion with strong evidence of full-length insertion and 
# subset info file
rownames(FullL1Info)[which(FullL1Info$Max3P >= 6020)]
idxFull2 <- which(rownames(FullL1Info) %in% FilesWithReads[idxFull])
FullL1InfoSubset <- FullL1Info[idxFull2,]

###############
# Align transduced sequences
###############

# Align transduced sequences and keep info about each insertion
AlignmentFilesTD5P <- c()
AlignmentFilesTD3P <- c()
FastaFiles_5P <- c()
FastaFiles_3P <- c()
ReadInfoListSubset        <- ReadInfoList[idxPotentialFull %in% idxFull2]
names(ReadInfoListSubset) <- FilesWithReads[idxFull]
for (i in idxFull2){
  
  # Get position of current entry in ReadListPerPeak
  x <- idx5P[i]
  
  # Get reads mapped to L1, retain only primary reads and get the number of bp
  # clipped on the left
  RL <- MatchedReadList[[which(i == idxPotentialFull)]]
  
  # Get Info and match reads to info
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  
  
  # Create a file name for output fasta file
  FastaFile5P <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",
                       FullL1Info$idx[i], "_TDS_5P.fas", sep = "")
  FastaFile3P <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",
                       FullL1Info$idx[i], "_TDS_3P.fas", sep = "")
  FastaFiles_5P <- c(FastaFiles_5P, FastaFile5P)
  FastaFiles_3P <- c(FastaFiles_3P, FastaFile3P)
  
  # Determine which sequence has enough clipped
  idxEnoughClipped <- which(
    pmax(L1Info$LeftClipped_G, L1Info$RightClipped_G) > MinClipAlign &
      L1Info$L1pos < 500)
  
  # Loop through reads and collect parts of reads that are transduced sequences
  TdS5PList   <- vector(mode = "list", length = length(idxEnoughClipped))
  TdS3PList   <- vector(mode = "list", length = length(idxEnoughClipped))
  idxL5P <- 0
  idxL3P <- 0
  for (y in idxEnoughClipped) {
    
    # Get sequence of current read and the number of left and right clipped bps
    Seq   <- RL$GenomeReadList$seq[y]
    if (L1Info$TD5Pstart[y] < L1Info$TD5Pend[y]){
      TdS5P   <- subseq(Seq, L1Info$TD5Pstart[y], L1Info$TD5Pend[y])
      idxL5P  <- idxL5P + 1
      TdS5PList[[idxL5P]] <- s2c(as.character(TdS5P))
    }
    if (L1Info$TD3Pstart[y] < L1Info$TD3Pend[y]){
      TdS3P   <- subseq(Seq, L1Info$TD3Pstart[y], L1Info$TD3Pend[y])
      idxL3P  <- idxL3P + 1
      TdS3PList[[idxL3P]] <- s2c(as.character(TdS3P))
    }
  }
  
  # Save transduced sequences as fasta file
  names(TdS5PList) <- RL$GenomeReadList$qname[idxEnoughClipped]
  names(TdS5PList) <- gsub("/", "_", names(TdS5PList))
  write.fasta(TdS5PList, names = names(TdS5PList), file.out = FastaFile5P)
  
  names(TdS3PList) <- RL$GenomeReadList$qname[idxEnoughClipped]
  names(TdS3PList) <- gsub("/", "_", names(TdS3PList))
  write.fasta(TdS3PList, names = names(TdS3PList), file.out = FastaFile3P)
  
  # Aligning transduced sequences as fasta file
  if (length(TdS5PList) > 1){
    AlignedFile5P <- gsub(".fas", "_aligned.fas", FastaFile5P)
    cat("Aligning", FastaFile5P, "\n")
    run_MUSCLE(InputPath = FastaFile5P, OutputPath = AlignedFile5P) 
    AlignmentFilesTD5P <- c(AlignmentFilesTD5P, AlignedFile5P)
  }
  if (length(TdS3PList) > 1){
    AlignedFile3P <- gsub(".fas", "_aligned.fas", FastaFile3P)
    cat("Aligning", FastaFile3P, "\n")
    run_MUSCLE(InputPath = FastaFile3P, OutputPath = AlignedFile3P) 
    AlignmentFilesTD3P <- c(AlignmentFilesTD3P, AlignedFile3P)
  }
}

# Loop through fasta files and create consensus sequences for transduced sequences
TDConsensus_5P <- lapply(FastaFiles_5P, function(x){
  AlignedFile <- gsub(".fas", "_aligned.fas", x)
  if(file.exists(AlignedFile)){
    Alignment <- read.dna(AlignedFile, format = "fasta", as.character = T, as.matrix = T)
    ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
      Nucs <- z[z != "-"]
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    })
  } else {
    ConsensusSeq <- read.dna(x, format = "fasta", as.character = T, as.matrix = T)
  }
  ConsensusSeq
})
names(TDConsensus_5P) <- FilesWithReads[idxFull]
TDConsensus_3P <- lapply(FastaFiles_3P, function(x){
  AlignedFile <- gsub(".fas", "_aligned.fas", x)
  if(file.exists(AlignedFile)){
    Alignment <- read.dna(AlignedFile, format = "fasta", as.character = T, as.matrix = T)
    ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
      Nucs <- z[z != "-"]
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    })
  } else {
    ConsensusSeq <- read.dna(x, format = "fasta", as.character = T, as.matrix = T)
  }
  ConsensusSeq
})
names(TDConsensus_3P) <- FilesWithReads[idxFull]

# Loop through potential full-length insertions and 
AlignmentFiles <- c()
i <- 68
for (i in idxFull2){
  print(i)
  # Index for transduced sequences
  idxTD <-  which(i == idxFull2)
  
  # Get Info and matched read list
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  
  RL     <- MatchedReadList[[which(i == idxPotentialFull)]]$GenomeReadList 
  
  # Create a file name for output fasta file
  FastaFile <- paste(FastaFilePath, FullL1Info$chromosome[i], "_",FullL1Info$idx[i],
                     ".fas", sep = "")
  
  # Get insertion location and create a reference sequence with sequence
  # flanking L1 and the consensus L1
  InsLoc  <- FullL1Info$L1InsertionPosition.median[i]
  L1Str   <- FullL1Info$L1Strand[i]
  L1start <- min(L1Info$L1pos)
  L1ConsSeqLocal <- L1CharV[L1start:length(L1CharV)]
  if (L1Str == "-") {
    L1ConsSeqLocal <- L1CharV_RV[1:(length(L1CharV_RV) - L1start + 1)]
    LeftTDs  <-  TDConsensus_3P[[idxTD]]
    RightTDs <-  TDConsensus_5P[[idxTD]]
  } else {
    LeftTDs  <-  TDConsensus_5P[[idxTD]]
    RightTDs <-  TDConsensus_3P[[idxTD]]
  }
  TransDL    <- FullL1Info$L15PTransdSeq.median[i]
  LeftFlankR <- GRanges(FullL1Info$chromosome[i], 
                        IRanges(start = InsLoc - FlankSize + 1,
                                end = InsLoc))
  RightFlankR <- GRanges(FullL1Info$chromosome[i], IRanges(start = InsLoc + 1,
                                                           end = InsLoc + FlankSize))
  LeftFlankSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, LeftFlankR)
  RightFlankSeq   <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RightFlankR)
  LeftFlankCharv  <- s2c(as.character(LeftFlankSeq))
  RightFlankCharv <- s2c(as.character(RightFlankSeq))
  
  RefSeqCharV     <- c(LeftFlankCharv, "N", LeftTDs, "N", 
                       L1ConsSeqLocal, "N", RightTDs, "N", RightFlankCharv)
  RefSeqCharV     <- toupper(RefSeqCharV)
  
  # Get reads mapped to the genome on current locus
  LRClippedG <- sapply(RL$cigar, NrClippedFromCigar)
  idxEnoughClipped <- which(
    pmax(LRClippedG[1,], LRClippedG[2,]) > MinClipAlign)
  
  # Loop through reads and collect parts of reads that should align with L1 and
  # its flanking sequence
  SeqList   <- vector(mode = "list", length = length(idxEnoughClipped) + 1)
  SeqList[[1]]     <- RefSeqCharV
  idxL <- 1
  for (y in idxEnoughClipped) {
    
    # Get sequence of current read and the number of left and right clipped bps
    Seq        <- RL$seq[y]
    LRClipped  <- LRClippedG[,y]
    
    # If more is clipped on the left take the sequence from the start into the
    # flank extending to the right of the cilpper part
    SeqStart <- switch(as.character(L1Info$Scenario[y]), 
                       '0' = max(1, L1Info$TD3Pstart[y] - FlankSize), 
                       '1' = max(1, L1Info$TD5Pstart[y] - FlankSize),
                       '2' = 1,
                       '3' = 1)
    SeqEnd <- switch(as.character(L1Info$Scenario[y]), 
                     '0' = width(Seq), 
                     '1' = width(Seq),
                     '2' = min(width(Seq), L1Info$TD5Pend[y] + FlankSize),
                     '3' = min(width(Seq), L1Info$TD3Pend[y] + FlankSize))
    Seq2Add  <- subseq(Seq[[1]], SeqStart, SeqEnd)
    idxL <- idxL + 1
    SeqList[[idxL]] <- s2c(as.character(Seq2Add))
  }
  names(SeqList) <- c("RefSeq", RL$qname[idxEnoughClipped])
  names(SeqList) <- gsub("/", "_", names(SeqList))
  write.fasta(SeqList, names = names(SeqList), file.out = FastaFile)
  
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
  
  # Read alignment and reorder so that reference sequence is first
  Alignment <- read.fasta(AlignedFile)
  idxRefSeq <- which(names(Alignment) == "RefSeq")
  AlignmentReordered <- Alignment
  AlignmentReordered[[1]] <- Alignment[[idxRefSeq]]
  names(AlignmentReordered)[1] <- names(Alignment)[idxRefSeq]
  AlignmentReordered[[idxRefSeq]] <- Alignment[[1]]
  names(AlignmentReordered)[idxRefSeq] <- names(Alignment)[1]
  write.fasta(AlignmentReordered, names = names(AlignmentReordered), file.out = AlignedFile)
  
  # Keep track of alignment file names
  AlignmentFiles <- c(AlignmentFiles, AlignedFile)
}

# Loop through alignment files and analyze regions of 
AlignFile <- AlignmentFiles[1]
AlignInfoList <- lapply(AlignmentFiles, function(AlignFile){
  
  # Extract plot title from file name
  Fsplit1 <- strsplit(AlignFile, "/")[[1]]
  Fsplit2 <- strsplit(Fsplit1[length(Fsplit1)], "_")[[1]]
  Title <- paste(Fsplit2[1:2], collapse = "-")
  
  # Read alignment and extract infor about alignment quality
  Alignment <- read.dna(AlignFile, format = "fasta", as.character = T, as.matrix = T)
  blnSame   <- sapply(1:ncol(Alignment), function(x){
    any(Alignment[2:nrow(Alignment),x] == Alignment[1,x])
  })
  blnSameTS <- as.ts(1*blnSame)
  SmoothedSame <- sma(blnSameTS, order = 20)
  nPos <- which(Alignment[1,] %in% c("n", "N"))
  
  # Extract consensus sequence
  Alignment <- Alignment[ , Alignment[1,] != "-"]
  # ConsensusSeq <- apply(Alignment, 2, FUN = function(z) {
  #   Nucs <- z[z != "-"]
  #   if (z[1] != "n"){
  #     NucFreq <- table(Nucs)
  #     names(NucFreq)[which.max(NucFreq)]
  #   } else {
  #     z[1]
  #   }
  # })
  ConsensusSeq <- sapply(1:ncol(Alignment), function(v) {
    z <- Alignment[,v]
    Nucs <- z[z != "-"]
    if (z[1] != "n" & v > min(nPos) & v < max(nPos)){
      NucFreq <- table(Nucs)
      names(NucFreq)[which.max(NucFreq)]
    } else {
      z[1]
    }
  })
  ConsensNPos <- which(ConsensusSeq %in% c("n", "N"))
  
  
  list(nPos = nPos, SmoothedSame = SmoothedSame, Title = Title,
       ConsensusSeq = ConsensusSeq, ConsensNPos = ConsensNPos)
})

# Get consensus sequence of flanking sequence and insertion
L1WithFlank <- lapply(AlignInfoList, function(x) x$ConsensusSeq)
names(L1WithFlank) <- sapply(AlignInfoList, function(x) x$Title)
names(L1WithFlank) <- gsub("-", "_", names(L1WithFlank))
write.fasta(L1WithFlank, names = names(L1WithFlank),
            file.out = L1InsertionFastaPath)

# Write reconstructed L1 insertions as individual fasta files so that reads
# can be aligned to them
for (i in 1:length(L1WithFlank)){
  FileNameL1 <- paste(FastaFilePath, "L1WithFlank", FlankSize, 
                      "bp_",names(L1WithFlank)[i], ".fas", sep = "")
  FileNameTdSeq <- paste(FastaFilePath, "TdSeq", FlankSize, 
                         "bp_",names(L1WithFlank)[i], ".bed", sep = "")
  FileNameL1pos <- paste(FastaFilePath, "L1pos", FlankSize, 
                         "bp_",names(L1WithFlank)[i], ".bed", sep = "")
  write.fasta(L1WithFlank[[i]], names = names(L1WithFlank)[i], FileNameL1)
  NPos      <- AlignInfoList[[i]]$ConsensNPos
  TransdPos <- data.frame(rep(names(L1WithFlank)[i], 2), NPos[c(1,3)], 
                          NPos[c(2,4)])
  L1Pos <- data.frame(names(L1WithFlank)[i], NPos[2], NPos[3])
  write.table(TransdPos, file = FileNameTdSeq, col.names = F, quote = F, 
              row.names = F, sep = "\t")
  write.table(L1Pos, file = FileNameL1pos, col.names = F, quote = F, 
              row.names = F, sep = "\t")
}

# Plot alignment quality
par(mfrow = c(3, 2), mar = c(3, 2, 3, 0.5), oma = c(2, 4, 1, 1))
for (i in 1:length(AlignInfoList)){
  SmoothedSame <- AlignInfoList[[i]]$SmoothedSame
  nPos <- AlignInfoList[[i]]$nPos
  Title <- AlignInfoList[[i]]$Title
  plot(SmoothedSame$fitted, type = "l", xlab = "", ylab = "", main = Title,
       ylim = c(0, 1.1))
  segments(x0 = nPos, y0 = 0, y1 = 1, col = "red", lwd = 2)
}
mtext("Proportion of bp same as in reference", 2, outer = T, line = 2)
CreateDisplayPdf("D:/L1polymORF/Figures/InsertionAlignment.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 10, width = 7)

# Get list of sequences with insertions
InsertionSeqs <- lapply(idxFull2, function(i){
  
  # Get Info
  L1Info <- ReadInfoList[[which(i == idxPotentialFull)]]  
  
  # Get insertion location and create a reference sequence with sequence
  # flanking L1 and the consensus L1
  InsLoc     <- FullL1Info$L1InsertionPosition.median[i]
  L1Str      <- NewStrand[i == idxFull2]
  L1ConsSeqLocal <- L1CharV
  if (L1Str == "-") L1ConsSeqLocal <- L1CharV_RV
  LeftFlankR <- GRanges(FullL1Info$chromosome[i], 
                        IRanges(start = InsLoc - FlankSize + 1,
                                end = InsLoc))
  RightFlankR <- GRanges(FullL1Info$chromosome[i], IRanges(start = InsLoc + 1,
                                                           end = InsLoc + FlankSize))
  LeftFlankSeq    <- getSeq(BSgenome.Hsapiens.UCSC.hg19, LeftFlankR)
  RightFlankSeq   <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RightFlankR)
  LeftFlankCharv  <- s2c(as.character(LeftFlankSeq))
  RightFlankCharv <- s2c(as.character(RightFlankSeq))
  c(LeftFlankCharv, L1ConsSeqLocal, RightFlankCharv)
})
SeqNames <- paste(FullL1InfoSubset$chromosome, FullL1InfoSubset$idx, sep = "_")
write.fasta(InsertionSeqs, names = SeqNames, file.out = "D:/L1polymORF/Data/L1InsertionSequences.fas")
InsertionSummaries <- FullL1InfoSubset[,c("chromosome", "L1InsertionPosition.median")]
colnames(InsertionSummaries)[colnames(InsertionSummaries) == "L1InsertionPosition.median"] <- "L1InsertionPosition"
InsertionSummaries$SeqName  <- SeqNames
InsertionSummaries$L1Strand <- NewStrand
InsertionSummaries$TransD5P <- c("Yes", "No", "No", "Yes",  "Yes")
InsertionSummaries$TransD3P <- c("No", "Yes", "No", "No",  "No")
write.csv(InsertionSummaries, file = "D:/L1polymORF/Data/L1InsertionSummary.csv",
          row.names = F)

# Create genomic ranges for L1 insertions and compare with known L1 ranges
GRL1Capture <- makeGRangesFromDataFrame(FullL1Info)
GRL1Capture100 <- resize(GRL1Capture, 100, fix = "center")
L1CatalogGR100 <- resize(L1CatalogGR, 100, fix = "center")
GRL1Ins1000G100 <- resize(GRL1Ins1000G, 100, fix = "center")
if (any(c(sum(overlapsAny(GRL1Capture100, L1CatalogGR100)),
          sum(overlapsAny(GRL1Capture100, GRL1Ins1000G100))) > 0)){
  cat("Some newly found L1 insertion overlap with catalog or 1000 Genome elements!\n")
}

###############################################
#                                             #
#       Assemble sequences of L1              #
#             insertions                      #
#                                             #
###############################################

# Get sequences of a particular insertion
NewL1Seq <- sapply(1:length(idx5P), function(i){
  x <- idx5P[i]
  cat("Processing insertion", i, "of", length(idx5P), "\n")
  RL <- ReadListPerPeak[[x]]
  idxReads <- which(RL$pos < FivePEnd)
  
  # Loop through reads and create a matrix of aligned reads
  ReadMat  <- sapply(idxReads, function(i) {
    SeqV   <- SeqFromCigar(RL$cigar[i], RL$seq[i])
    SeqEnd  <- min(FivePEnd - RL$pos[i] + 1, length(SeqV))
    Prepend <- rep("-", RL$pos[i] - 1)
    NrAppend <- FivePEnd - RL$pos[i] + 1 - SeqEnd
    Append  <- rep("-", NrAppend)
    SV <- c(Prepend, SeqV[1:SeqEnd], Append)
  })
  
  # Create a consensus sequence
  ConsensSeq  <- L1ConsSeq[1:FivePEnd]
  ConsensProp <- rep(1, FivePEnd)
  AllDel <- sapply(1:nrow(ReadMat), function(x) all(ReadMat[x,] == "-"))
  ConsensSeq[AllDel] <- "-"
  for(x in which(!AllDel)){
    NucCount    <- table(ReadMat[x,])
    NucCount    <- NucCount[names(NucCount) != "-"]
    ConsensSeq[x]  <- names(NucCount)[which.max(NucCount)]
    ConsensProp[x] <- max(NucCount)/sum(NucCount)
  }
  
  # Replace nucleotides that are variable 
  idxReplace <- which(ConsensProp < MinPolyProp)
  ConsensSeq[idxReplace] <- L1ConsSeq[idxReplace]
  ConsensSeq
})
colnames(NewL1Seq) <- paste(FullL1Info$chromosome, FullL1Info$start,
                            FullL1Info$end, "New", sep = "_")

# Calculate the number of nucleotide diefferences between different L1 stumps
DiffMat <- sapply(1:ncol(NewL1Seq), function(x) {
  cat("Calculating differences for sequence", x, "of", ncol(NewL1Seq), "\n")
  sapply(1:ncol(NewL1Seq), function(y){
    sum(NewL1Seq[ , x] != NewL1Seq[ , y])
  })
})

# Determine indices of sequences that are identical with another
diag(DiffMat) <- NA
max(DiffMat, na.rm = T)
min(DiffMat, na.rm = T)
mean(DiffMat, na.rm = T)
sum(DiffMat == 0, na.rm = T)
idxIdentical <- which(DiffMat == 0, arr.ind = T)
idxIdentical <- unique(as.vector(idxIdentical))

# Count differences to consensus sequence
Diff2L1Consens <- sapply(1:ncol(NewL1Seq), function(x) {
  sum(NewL1Seq[ , x] != L1ConsSeq[1:FivePEnd])
})
hist(Diff2L1Consens, breaks = seq(-5, 1005, 5))
sum(Diff2L1Consens == 0)
mean(Diff2L1Consens)
NrDel <- sapply(1:ncol(NewL1Seq), function(x) {
  sum(NewL1Seq[ , x] == "-")
})
rownames(FullL1Info)[which.max(NrDel)]

# Get L1 5' sequences that are non-reference
L1Catalogue <- LiftOverList$L1CatalogWithHG19
blnNonRef   <- (L1Catalogue$end_HG19 - L1Catalogue$start_HG19) <= 5000
idxNonRef   <- which(blnNonRef)
L1SeqNonRef <- L1Catalogue$L1Seq[idxNonRef]

# Turn all L1 seq 
# Pattern <- paste(L1ConsSeq[100:200], collapse = "")
# pMatch <- vmatchPattern(Pattern, L1SeqNonRef, max.mismatch = 5)
# sapply(pMatch, length)
L1NonRef <- sapply(L1SeqNonRef, function(x) strsplit(x, "")[[1]][1:FivePEnd])
colnames(L1NonRef) <- paste(L1Catalogue$Chromosome[idxNonRef], 
                            L1Catalogue$start_HG19[idxNonRef],
                            L1Catalogue$end_HG19[idxNonRef], 
                            L1Catalogue$Accession[idxNonRef], sep = "_")


# Put all the different L1 sequences together, remove consistent deletions and
# save as fasta file
L1StumpRef <- cbind(NewL1Seq, L1NonRef, L1RefSeqMat)
L1StumpList <- lapply(1:ncol(L1StumpRef), function(x) {
  c(L1StumpRef[L1StumpRef[,x] != "-",x], L1ConsSeq[(FivePEnd + 1):length(L1ConsSeq)])
})
cat("*****   Saving new L1 references as file", NewL1RefOutPath, "  *****\n")
write.fasta(L1StumpList, colnames(L1StumpRef), NewL1RefOutPath)

# Save info about insertions
cat("*****   Saving image to", ResultPath, "  *****\n")
save.image(ResultPath)
