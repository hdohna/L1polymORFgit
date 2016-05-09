##############################################
#
# General description:
#
#   The following script reads in a catalogue of functional L1s and uses bwa to
#   find flanking regions of L1 that can only be found on BC clones/fosmids

# Input:
#
#    D:/L1polymORF/Data/L1catalogue.csv: table with functional L1HS
#   

# Output:
#   
#    L1HS_repeat_table.csv: csv file with all L1HS repeats
#    L1Sequences_reptab.fas: fasta file with all L1 sequences

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load packages
library(Rsamtools)
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(seqinr)

# Specify data path
DataPath <- "/home/hzudohna/L1polymORF/Data/"
#DataPath <- "D:/L1polymORF/Data/"

# Specify commands
IndexCommand <- '/home/txw/bwa/bwa-0.7.12/bwa index'
AlignCommand <- '/home/txw/bwa/bwa-0.7.12/bwa mem'

# Indicators for different analysis steps
blnWriteFasta <- F
blnBuildIndex <- F
blnConvertSam2Bam <- T

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Read catalogue
L1CataloguePath <- paste(DataPath, "L1Catalogue_Sat_May_07_15-15-31_2016.csv", sep = "")
L1Catalogue     <- read.csv(L1CataloguePath, as.is = T)

#######################################
#                                     #
#     Write reference seqs            #
#                                     #
#######################################

# Get indices of entries whose flanks need to be searched
idxFlanksToBeSearched <- which(is.na(L1Catalogue$start_HG38) & 
                                 !is.na(L1Catalogue$L1SeqFlank5p))

# Get chromosomes that need to be searched and write a reference sequence
# for each
ChrsToBeSearched <- unique(L1Catalogue$Chromosome[idxFlanksToBeSearched])
ChrPathVect <- c()
cat("******  Writing out a reference sequence per chromosome  ***********\n")
for (Chr in ChrsToBeSearched){
  cat("Writing chromosome", Chr, "\n")
  ChrPath <- paste(DataPath, Chr, ".fas", sep = "")
  if (blnWriteFasta){
    ChrSeq  <- BSgenome.Hsapiens.UCSC.hg38[[Chr]]
    write.fasta(s2c(as.character(ChrSeq)), names = Chr, file.out = ChrPath)
  }
  ChrPathVect <- c(ChrPathVect, ChrPath)
}
names(ChrPathVect) <- ChrsToBeSearched

# Run command to create indices for each chromosome
if(blnBuildIndex){
  for (ChrPath in ChrPathVect){
    CmdIndex <- paste(IndexCommand, ChrPath)
    system(CmdIndex)
  }
  
}

# Write each flanking sequence as fastq file and run an alignment command
cat("******  Writing out fastq file per L1 and aligning them ***********\n")
for (i in idxFlanksToBeSearched){
  cat("\n\n Writing fastq file", i, "\n")
  Chr <- L1Catalogue$Chromosome[i]
  ReadList <- list(seq  = c(L1Catalogue$L1SeqFlank5p[i],
                            L1Catalogue$L1SeqFlank3p[i]),
                   qual = c(paste(rep("~", nchar(L1Catalogue$L1SeqFlank5p[i])),
                                  collapse = ""),
                            paste(rep("~", nchar(L1Catalogue$L1SeqFlank5p[i])),
                                           collapse = "")),
                   rname = c(Chr, Chr),
                   qname = c(Chr, Chr),
                   pos = c("",""))
  FilePrefix <- paste("L1Flank", L1Catalogue$Accession[i], sep = "_")
  WriteFastqAndSample(ReadList, DataPath, FilePrefix = FilePrefix)
  FastqPath <- paste(DataPath, FilePrefix, ".fastq", sep = "")
  OutFile <- paste(DataPath, FilePrefix, ".sam", sep = "")
  CmdLine <- paste(AlignCommand,  ChrPathVect[Chr], FastqPath)
  CmdLine <- paste(CmdLine, OutFile, sep = " > ")
  system(CmdLine)
}

# Get paths to all sam files
FilePaths <- list.files(DataPath, pattern = "L1Flank", full.names = T)
SamFilePaths <- grep(".sam", FilePaths, value = T)
for (SamFile in SamFilePaths){
  
  # Get accession number from file name
  AccNr  <- strsplit(SamFile, "_")[[1]][2]
  AccNr  <- strsplit(AccNr, "\\.")[[1]][1]
  
  # Determine matching row in catalogue, chromsome and create GRanges object
  idxRow <- which(L1Catalogue$Accession == AccNr)
  Chrom  <- L1Catalogue$Chromosome[idxRow]
  GR     <- GRanges(seqnames = Chrom, IRanges(start = 1, 
               end = length(BSgenome.Hsapiens.UCSC.hg38[[Chrom]])))
  BamFile <- gsub(".sam", ".bam", SamFile)
  if (blnConvertSam2Bam & file.exists(BamFile)){
    BamFilePrefix <- substr(BamFile, 1, nchar(BamFile) - 4)
    asBam(SamFile, BamFilePrefix, overwrite = T)
    sortBam(BamFile, BamFile)
    indexBam(BamFile)
    
  }
  if (file.exists(BamFile)){
    cat("Extracting reads from", BamFile, "\n")
    ReadRanges <- extractReads(BamFile, GR)
    NewStart   <- min(end(ReadRanges))
    NewEnd     <- max(start(ReadRanges))
    length(ReadRanges)
    if (length(ReadRanges) == 2 & ((NewEnd - NewStart) < 6100) &
        is.na(L1Catalogue$end_HG38[idxRow])){
      cat("Adding start and end for", AccNr, "\n")
      L1Catalogue$start_HG38[idxRow] <- NewStart
      L1Catalogue$end_HG38[idxRow] <- NewEnd
    }
  } else {
    cat("Could not find bam file", BamFile, "\n")
  }
}

# Write updated catalogue out
SplitPath <-  strsplit(L1CataloguePath, "_")[[1]]
CataloguePath <- paste(c(SplitPath[1], "Updated", SplitPath[-1]), collapse = "_")
cat("Saving new L1 catalogue file", CataloguePath, "\n")
write.csv(L1Catalogue, CataloguePath)
