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
library(BSgenome.Hsapiens.UCSC.hg38)
#library(seqinr)

# Specify data path
DataPath <- "/home/hzudohna/L1polymORF/Data/"
#DataPath <- "D:/L1polymORF/Data/"

# Specify commands
IndexCommand <- '/home/txw/bwa/bwa-0.7.12/bwa index'
AlignCommand <- '/home/txw/bwa/bwa-0.7.12/bwa mem'

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Read catalogue
L1CataloguePath <- paste(DataPath, "L1Catalogue.csv", sep = "")
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
  ChrSeq  <- BSgenome.Hsapiens.UCSC.hg38[[Chr]]
  ChrPath <- paste(DataPath, Chr, ".fas", sep = "")
#  write.fasta(s2c(as.character(ChrSeq)), names = Chr, file.out = ChrPath)
  ChrPathVect <- c(ChrPathVect, ChrPath)
}
names(ChrPathVect) <- ChrsToBeSearched

# Run command to create indices for each chromosome
for (ChrPath in ChrPathVect){
  CmdIndex <- paste(IndexCommand, ChrPath)
  system(CmdIndex)
}

# Write each flanking sequence as fastq file and run an alignment command
cat("******  Writing out fastq file per L1 and aligning them ***********\n")
for (i in idxFlanksToBeSearched){
  cat("Writing fastq file", i, "\n")
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
