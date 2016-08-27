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
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(Rsamtools)
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg38)
library(seqinr)

# Specify data path and file name for L1 catalog
DataPath <- "/home/hzudohna/L1polymORFgit/Data/"
L1CatalogFile <- "L1Catalog_Wed_Aug_10_17-32-20_2016.csv"

# Specify commands
IndexCommand <- c('module load bwa', 'bwa index')
AlignCommand <- c('module load bwa', 'bwa mem')
SamSortCommand <- c('module load samtools', 
   'InFile=$1', 
   'TmpFile=${InFile%.sam}.tmp', 
   'SortedBamFile=${InFile%.sam}.sorted.bam', 
   'samtools view $InFile -b -h | samtools sort -o $SortedBamFile -T $TmpFile -'
)

# Indicators for different analysis steps
blnWriteFasta <- F
blnBuildIndex <- F
blnAlnFlanks  <- F
blnConvertSam2Bam <- F
blnIndexBam <- F

# Define function to switch strand
StrandSwitch <- function(Strand) switch(Strand, '+' = "-", '-' = "+")

# Define standard header lines
qsubHeaderLines = c('#! /bin/sh', 
                    '#$ -N TEST', 
                    '#$ -cwd', '#', 
                    '#$ -l h_rt=1:00:00', 
                    '#$ -j y', '#$ -S /bin/bash', '#', '')

#######################################
#                                     #
#     Read data                       #
#                                     #
#######################################

# Load chromosome lengths
load('/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/ChromLengthsHg38.Rdata')
ChromLengths <- ChromLengthsHg38

# Read catalogue
L1CataloguePath <- paste(DataPath, L1CatalogFile, sep = "")
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
    CmdIndex <- c(IndexCommand[1], paste(IndexCommand[2], ChrPath))
    Chr <- strsplit(ChrPath, "/")[[1]]
    Chr <- Chr[length(Chr)]
    ScriptFile <- paste('/home/hzudohna/qsub_index', Chr, sep = "_")
    scriptName <- paste('indextmp', Chr, sep = "_")
    CreateAndCallqsubScript(file = ScriptFile, 
       qsubHeaderLines = qsubHeaderLines, qsubCommandLines = CmdIndex,
       scriptName = scriptName)
  }
}

# Write each flanking sequence as fastq file and run an alignment command
if (blnAlnFlanks){
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
    CmdLine <- paste(AlignCommand[2],  ChrPathVect[Chr], FastqPath)
    CmdLine <- paste(CmdLine, OutFile, sep = " > ")
    CmdLine <- c(AlignCommand[1], CmdLine)
    ScriptFile <- paste('/home/hzudohna/qsub_align', L1Catalogue$Accession[i], sep = "_")
    scriptName <- paste('bwa', L1Catalogue$Accession[i], sep = "_")
    CreateAndCallqsubScript(file = ScriptFile, 
                            qsubHeaderLines = qsubHeaderLines,
                            qsubCommandLines = CmdLine,
                            scriptName = scriptName)
  }
}

# Get paths to all sam files
FilePaths <- list.files(DataPath, pattern = "L1Flank", full.names = T)
SamFilePaths <- grep(".sam", FilePaths, value = T)
for (SamFile in SamFilePaths){
  if (blnConvertSam2Bam & file.exists(SamFile)){
    cat("Conerting sam file", SamFile, "\n")
    # Get accession number from file name
    AccNr  <- strsplit(SamFile, "_")[[1]][2]
    AccNr  <- strsplit(AccNr, "\\.")[[1]][1]
    
    ScriptFile <- paste('/home/hzudohna/qsub_SamSort', AccNr, sep = "_")
    scriptName <- paste('SamSorttmp', AccNr, sep = "_")
    CreateAndCallqsubScript(file = ScriptFile, 
       qsubHeaderLines = qsubHeaderLines, 
       qsubCommandLines = SamSortCommand,
       scriptName = scriptName, Args = SamFile)
  }
}

# Create vector of bam files
BamFilePaths <- gsub(".sam", ".sorted.bam", SamFilePaths)

# Index each bam file
if (blnIndexBam){
  for (BamFile in BamFilePaths){
    if (file.exists(BamFile)){
      cat("Indexing bam file", BamFile, "\n")
      CmdLine <- paste("qsub /home/hzudohna/pbs_SamIndex", BamFile)
      system(CmdLine)
    }
  }
}

# Waiting loop for indexing to be done
# QStat <- system("qstat", intern = T)
# jobCount <- length(grep("SamIndex", QStat))
# while (jobCount > 0){
#   cat("Waiting for index jobs to finish ...\n")
#   QStat <- system("qstat", intern = T)
#   jobCount <- length(grep("SamIndex", QStat))
#   Sys.sleep(30)
# }


for (BamFile in BamFilePaths){
  if (file.exists(BamFile)){
    cat("Extracting reads from", BamFile, "\n")
    # Get accession number from file name
    AccNr  <- strsplit(BamFile, "_")[[1]][2]
    AccNr  <- strsplit(AccNr, "\\.")[[1]][1]

    # Determine matching row in catalogue, chromsome and create GRanges object
    idxRow <- which(L1Catalogue$Accession == AccNr)
    Chrom  <- L1Catalogue$Chromosome[idxRow]
    GR     <- GRanges(seqnames = Chrom, IRanges(start = 1,
                 end = ChromLengths[Chrom]))
    ReadRanges <- extractReads(BamFile, GR, param = readParam(forward = T))
    if (length(ReadRanges) == 0){
      ReadRanges <- extractReads(BamFile, GR)
      if (length(ReadRanges) > 0){
        cat("Correcting strand for", AccNr, "\n")
        L1Catalogue$strand_ClonetoRef[idxRow] <- '-'
      }
    } else {
      L1Catalogue$strand_ClonetoRef[idxRow] <- '+'
    }   
    NewStart   <- min(end(ReadRanges))
    NewEnd     <- max(start(ReadRanges))
    length(ReadRanges)
    if (length(ReadRanges) == 2 & ((NewEnd - NewStart) < 6100) &
        is.na(L1Catalogue$end_HG38[idxRow])){
      cat("Adding start and end for", AccNr, "\n")
      L1Catalogue$start_HG38[idxRow] <- NewStart
      L1Catalogue$end_HG38[idxRow]   <- NewEnd
    }
  } else {
    cat("Could not find bam file", BamFile, "\n")
  }
}

# Update strandedness with respect to reference
L1Catalogue$strand_L1toRef <- L1Catalogue$Strand
blnStrandDiff <- L1Catalogue$strand_ClonetoRef == '-'
SwitchedStrands <- sapply(L1Catalogue$Strand[which(blnStrandDiff)], 
                          function(x) StrandSwitch(x))
L1Catalogue$strand_L1toRef[which(blnStrandDiff)] <- 
  SwitchedStrands
L1Catalogue$strand_L1toRef[is.na(L1Catalogue$start_HG38)] <- NA

# Write updated catalogue out
SplitPath <-  strsplit(L1CataloguePath, "_")[[1]]
CataloguePath <- paste(c(SplitPath[1], "Updated", SplitPath[-1]), collapse = "_")
cat("Saving new L1 catalogue file", CataloguePath, "\n")
write.csv(L1Catalogue, CataloguePath)
