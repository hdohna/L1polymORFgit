##############################################
#
# General description:
#
#   The following script reads a L1 catalog and creates a separate alignment for
#   all reads that overlap the junction.

# Input:
#
#     BamFile: path to file that contains mapped reads
#     L1TableFileName: path to file that contains L1HS ranges in a table

# Output:
#   
#    : ...

##############################################

######                                      
# Source packages and set parameters  
######

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(ape)
library(ShortRead)
library(rtracklayer)
library(Rsamtools)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Files and folders
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
L1CataloguePath <- "D:/L1polymORF/Data/L1Catalog_Updated_Wed_Aug_10_17-32-20_2016.csv"
FastaFolder     <- "D:/L1polymORF/Data/"

ChainFile       <- 'D:/L1polymORF/Data/hg38ToHg19.over.chain'
BamFilePath     <- 'D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19withL1.sorted.bam'

# Specify the minimum read depth to create alignment
MinReadDepth <- 5
MinReadNr    <- 5
MinMapQ      <- 30

# Specify the offset from the insertion junction
JunctionOffset <- 50

#######                       
# Get L1 ranges                    
#######                                     

cat("Getting reference L1 ranges \n")

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1),]
sum(blnL1Mapped)
sum(blnAllele1, na.rm = T)
sum(blnInRef, na.rm = T)

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
    ChainFilePath = ChainFile)
L1CatalogGR <- LiftOverList$GRCatalogue_hg19
L1CatalogL1Mapped <- L1CatalogL1Mapped[LiftOverList$idxUniqueMapped, ]
width(L1CatalogGR)

# Read in L1 ranges (created in script CreateHg19NonRefL1BedRanges)
L1RangesDF   <- read.table("D:/L1polymORF/Data/L1Ranges_hg19_withNonRefL1.txt", as.is = T,
                           header = T)

# Turn L1 ranges into a GRanges object
L1Ranges <- makeGRangesFromDataFrame(L1RangesDF)
L1RangesLeft  <- flank(L1Ranges, width = 1,  ignore.strand = T)
L1RangesRight <- flank(L1Ranges, width = 1, start = F,
                       ignore.strand = T)
L1RangesLeft  <- shift(L1RangesLeft, -JunctionOffset)
L1RangesRight <- shift(L1RangesRight, JunctionOffset)

ReadRangeListLeft <- lapply(1:length(L1RangesLeft), function(x){
  extractReads(BamFilePath, L1RangesLeft[x], param = readParam(minq = MinMapQ))
})
ReadRangeListRight <- lapply(1:length(L1RangesRight), function(x){
  extractReads(BamFilePath, L1RangesRight[x], param = readParam(minq = MinMapQ))
})

# StuffToScan <- c("pos", "seq", "strand", "rname", "cigar")
# paramScanLeft  <- ScanBamParam(which = L1RangesLeft, what = StuffToScan)
# paramScanRight <- ScanBamParam(which = L1RangesRight, what = StuffToScan)
# ReadListLeft   <- scanBam(BamFilePath, param = paramScanLeft)
# ReadListRight  <- scanBam(BamFilePath, param = paramScanRight)

# # Match 
# L1RangeNames <- as.vector(seqnames(L1Ranges))
# blnRefRanges <- nchar(L1RangeNames) <= 5
# idxRefRanges <- which(blnRefRanges)
# L1RefRanges  <- L1Ranges[blnRefRanges]
# idxOverlaps  <- findOverlaps(L1RefRanges, L1CatalogGR, minoverlap = 6000)
# RangeAccession <- rep(NA, length(L1Ranges))
# RangeAccession[idxRefRanges[idxOverlaps@queryHits]] <- L1CatalogL1Mapped$Accession[idxOverlaps@subjectHits]
# RangeAccession[!blnRefRanges] <- L1RangeNames[!blnRefRanges]

#####                                   
# Analyze insertion flanks                       
#####

# # Initialize a result table
# ResultTable <- L1CatalogL1Mapped[,c("Accession", "Chromosome", 
#         "strand_L1toRef", "start_HG38", "end_HG38")]
# ResultTable$start_HG19  <- start(L1CatalogGR)
# ResultTable$end_HG19    <- end(L1CatalogGR)
# ResultTable$TotalNrReads  <- NA
# ResultTable$Coverage5P  <- NA
# ResultTable$Coverage3P  <- NA
# ResultTable$ReadWidth5P <- NA
# ResultTable$ReadWidth3P <- NA
# ResultTable$RelJunctPos5P <- NA
# ResultTable$RelJunctPos3P <- NA

ResultTable <- L1RangesDF
ResultTable$TotalNrReads  <- NA
ResultTable$CoverageLeft  <- sapply(ReadRangeListLeft, length)
ResultTable$CoverageRight <- sapply(ReadRangeListRight, length)
names(ResultTable)[names(ResultTable) == "name"] <- "Accession"
ResultTable$MaxCover <- pmax(ResultTable$CoverageLeft, ResultTable$CoverageRight)

# Write out results table
write.csv(ResultTable, "D:/L1polymORF/Data/L1PacBioResults.csv", row.names = F)


# # Create a vector for start and end of a range
# StartEnd  <- c("start", "end")
# Rshift    <- c(-JunctionOffset, JunctionOffset)
# 
# # Loop through L1s and create local alignments of all reads that intersect
# #  with the flank of a particular L1
# FastaFileNames <- c()
# for (i in 1:nrow(ResultTable)){
#   Acc <- ResultTable$Accession[i]
#   Str <- ResultTable$strand_L1toRef[i]
#   cat("Analyzing", Acc, "\n")
# 
#   # Get range to intersect with reads (shifted to the outside of L1)
#   FlankPos <- c(start(L1CatalogGR[i]), end(L1CatalogGR[i])) + Rshift 
#   R <- GRanges(seqnames = rep(L1CatalogL1Mapped$Chromosome[i], 2),
#                IRanges(start = FlankPos, end = FlankPos))
#   paramScan  <- ScanBamParam(which = R, 
#                   what = c("pos", "seq", "strand", "rname", "cigar"))
#   Reads      <- scanBam(BamFilePath, param = paramScan)
# 
#   NrReads <- sapply(Reads, function(x) length(x$seq))
#   Reads   <- Reads[NrReads > 0]
#   if (any(NrReads > 0)){
#     
#     # Retain only reads that actually intersect (bug in scanBam?)
#     Reads <- subsetReadListByOverlap(Reads, R)
#     
#     # Get sequences in common range 
#     Seqs      <- Reads[[1]]$seq
#     Strands   <- Reads[[1]]$strand
#     Positions <- Reads[[1]]$pos
#     Cigars    <- Reads[[1]]$cigar
#     if (length(Reads) > 1){
#         Seqs      <- c(Seqs, Reads[[2]]$seq)
#         Strands   <- c(Strands, Reads[[2]]$strand)
#         Positions <- c(Positions, Reads[[2]]$pos)
#         Cigars    <- c(Cigars, Reads[[2]]$cigar)
#     }
#     
#     # Create a unique ID and remove duplicates
#     UniqueID <- paste(Strands, Positions, Cigars, sep = "_")
#     blnDupl  <- duplicated(UniqueID)
#     Seqs     <- Seqs[!blnDupl]
#     Strands  <- Strands[!blnDupl]
# #    Seqs[Strands == 2] <- reverseComplement(Seqs[Strands == 2])
#     ResultTable$TotalNrReads[i] <- length(Seqs)
#     # Get the reference sequence in common range
#     GRL1 <- L1CatalogGR[i]
#     strand(GRL1) <- "+"
#     RefSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRL1)
#     RefSeq      <- as(RefSeq, "DNAStringSet")
#     
#     # Create a file name for fasta file
#     SeqCounter <- 1
#     if (length(Seqs) > MinReadNr){
#       for (idxCurrentSeq in 1:length(Seqs)){
#         SeqsWithRef <- c(RefSeq, Seqs[idxCurrentSeq])
#         FastaFileName <- paste(FastaFolder, "PacBioSeqs_", Acc, "_", SeqCounter, ".fas", 
#                                sep = "")
#         cat("Writing sequences to", FastaFileName, "\n")
#         writeFasta(SeqsWithRef, FastaFileName)
#         FastaFileNames <- c(FastaFileNames, FastaFileName)
#         SeqCounter <- SeqCounter + 1
#       }
#       
#     }
#     
#   }
# }


# Loop through fasta files and align them
for (FastaFile in FastaFileNames){
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
}

# Create vector of all alignment file names 
AlignmentFiles <- gsub(".fas", "_aligned.fas", FastaFileNames)

# Loop through L1s and create local alignments of single reads with reference
# sequence
# FastaFileNames <- c()
# for (i in 1:nrow(ResultTable)){
#   Acc <- ResultTable$Accession[i]
#   Str <- ResultTable$strand_L1toRef[i]
#   cat("Analyzing", Acc, "\n")
#   for (j in 1:length(StartEnd)){
#     if ((Str == "+" & j == 1) |  (Str == "-" & j == 2)){
#       L1Side <- '5P'
#     } else {
#       L1Side <- '3P'
#     }
#     
#     # Get range to intersect with reads (shifted to the out side of L1)
#     R <- resize(L1CatalogGR[i], 1, fix = StartEnd[j], ignore.strand = T)
#     strand(R) <- "*"
#     R <- shift(R, Rshift[j])
#     paramScan  <- ScanBamParam(which = R, what = scanBamWhat())
#     Reads      <- scanBam(BamFilePath, param = paramScan)
#     
#     # Retain only reads that actually intersect (bug in scanBam?)
#     if (length(Reads[[1]]$pos) > 0){
#       ReadGR <- GRanges(seqnames = Reads[[1]]$rname[1],
#                         ranges = IRanges(Reads[[1]]$pos, width = width(Reads[[1]]$seq)))
#       blnOverlap <- overlapsAny(ReadGR, R)
#       Reads[[1]] <- lapply(Reads[[1]], function(x) x[blnOverlap])
#     }
#     
#     # Update coverage at flank
#     NrReadCol <- paste('Coverage', L1Side, sep = "")
#     ResultTable[i, NrReadCol] <- length(Reads[[1]]$pos)
#     
#     cat("Start", start(R), "\n")
#     cat("Strand", Str, "\n")
#     cat("The coverage at", i, NrReadCol, "equals", length(Reads[[1]]$pos), "\n")
#     if (length(Reads[[1]]$pos) > MinReadDepth){
#       
#       # Get common start and end for all reads
#       StartAll <- min(Reads[[1]]$pos)
#       EndAll   <- max(Reads[[1]]$pos + width(Reads[[1]]$seq))
#       RelStart <- StartAll - Reads[[1]]$pos + 1
#       RelEnd   <- RelStart + EndAll - StartAll - 1
#       
#       # Fill in entry for position of junction within alignment
#       RelJunctPosCol <- paste("RelJunctPos", L1Side, sep = "")
#       ResultTable[i, RelJunctPosCol] <- start(R) - StartAll - Rshift[j]
#       
#       # Get sequences in common range 
#       Seqs <- Reads[[1]]$seq
#       cat("Getting sequences\n")
#       
#       # Get the reference sequence in common range
#       RefSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, 
#                        names = Reads[[1]]$rname[1], 
#                        start = StartAll, end = EndAll)
#       RefSeq      <- as(RefSeq, "DNAStringSet")
#       SeqsWithRef <- c(RefSeq, Seqs)
#       names(SeqsWithRef) <- c("Reference", Reads[[1]]$qname)
#       # Create a file name for fasta file
#       FastaFileName <- paste(FastaFolder, "PacBioSeqs_", Acc, "_", L1Side,
#                              "_", k, ".fas", sep = "")
#       cat("Writing sequences to", FastaFileName, "\n")
#       writeFasta(SeqsWithRef, FastaFileName)
#       FastaFileNames <- c(FastaFileNames, FastaFileName)
#       
      
      # for (k in 1:length(Seqs)){
      #   cat("Sequence", k, "with strand", Reads[[1]]$strand[k], "\n")
      #   #if (i == 20 & k == 7 & L1Side == '5P') browser()
      #   if (Reads[[1]]$strand[k] == 2){
      #     
      #     Seqs[k] <- reverseComplement(Seqs[k])
      #   }
      #   Seqs[k] <- as(Seqs[k][[1]][RelStart[k]:RelEnd[k]], "DNAStringSet")
      #   
      #   # Append reference sequence to current sequence
      #   SeqsWithRef <- c(RefSeq, Seqs[k])
      #   
      #   # Create a file name for fasta file
      #   FastaFileName <- paste(FastaFolder, "FlankSeqs_", Acc, "_", L1Side, 
      #                          "_", k, ".fas", sep = "")
      #   cat("Writing sequences to", FastaFileName, "\n")
      #   writeFasta(SeqsWithRef, FastaFileName)
      #   FastaFileNames <- c(FastaFileNames, FastaFileName)
      # }
      
      # Fill in entry for alignment width
  #     ReadWidthCol <- paste("ReadWidth", L1Side, sep = "")
  #     ResultTable[i, ReadWidthCol] <- width(SeqsWithRef)[1]
  #     
  #   }
  # }
# }


# Create character string to identify flanks
# FlankIDs <- sapply(AlignmentFiles, function(x) {
#   paste(strsplit(x, "_")[[1]][2:3], collapse = "_")
# })
FlankIDs <- sapply(AlignmentFiles, function(x) {
  strsplit(x, "_")[[1]][2]
})
FlankIDs <- unique(FlankIDs)

# Loop through flank IDs, collect all alignments per ID and append them to one
# alingment
NewAlignmentFiles <- c()
for (FlankID in FlankIDs){
  IDParts   <- strsplit(FlankID, "_")[[1]]
  ResultRow <- which(ResultTable$Accession == IDParts[1])
  RelJunctPosCol <- paste("RelJunctPos", IDParts[2], sep = "")
  JunctPos <- ResultTable[ResultRow, RelJunctPosCol]
  CurrentAlignments <- grep(FlankID, AlignmentFiles, value = T)
  LocalAlignment <- read.dna(CurrentAlignments[1], as.matrix = T, 
                             format = "fasta", as.character = T)
  RefSeq <- LocalAlignment[1, LocalAlignment[1, ] != "-", drop = F]
  NewAlignment <- t(sapply(CurrentAlignments, function(x){
    LocalAlignment <- read.dna(x, as.matrix = T, format = "fasta", as.character = T)
    LocalAlignment[2, LocalAlignment[1, ] != "-"]
  }))
  NewAlignment <- rbind(RefSeq, NewAlignment)
  OutputName  <- paste(FastaFolder, "PacBioCombinedAln_", FlankID, ".fas", sep = "")
  write.dna(NewAlignment, OutputName, format = "fasta")
  NewAlignmentFiles <- c(NewAlignmentFiles, OutputName)
  
  # NewAlignmentTrunc <- NewAlignment[,(JunctPos - 25):(JunctPos + 25)]
  # OutputNameTrunc  <- paste(FastaFolder, "FlankCombAlnTrunc_", FlankID, ".fas", sep = "")
  # write.dna(NewAlignmentTrunc, OutputNameTrunc, format = "fasta")
  
}

# Create consensus sequence
NewAlignmentFile <- NewAlignmentFiles[1]
NewAlignment     <- read.dna(NewAlignmentFile, format = "fasta", as.character = T, as.matrix = T)

L1Consensus <- NewAlignment[1,]
OtherSeqs   <- NewAlignment[2:nrow(NewAlignment),]
NoBlank     <- OtherSeqs != "-"
blnMajorityDifferent <- sapply(1:ncol(NewAlignment), function(x) {
  sum(OtherSeqs[NoBlank[,x],x] != L1Consensus[x]) / sum(NoBlank[,x]) > 0.8
})
NucFreqDiff <- lapply(which(blnMajorityDifferent), function(x) table(OtherSeqs[NoBlank[,x],x]))

NewNucs <- sapply(which(blnMajorityDifferent), function(i){
  NucFreqs <- table(OtherSeqs[NoBlank[,i],i])
  names(NucFreqs)[which.max(NucFreqs)]
})
L1PolymORFs <- data.frame(Position = which(blnMajorityDifferent), RefNuc = L1Consensus[blnMajorityDifferent], 
      NewNuc = NewNucs)
write.csv(L1PolymORFs, "D:/L1polymORF/Data/L1polymORFs_AC018980.csv", row.names = F)
write.fasta(L1Consensus, names = "AC018980", file.out = "D:/L1polymORF/Data/AssembledL1_AC018980.fas")
  
# boxplot(ResultTable$Coverage5P, ResultTable$Coverage3P)
# mean(ResultTable$Coverage5P)
# mean(ResultTable$Coverage3P)
