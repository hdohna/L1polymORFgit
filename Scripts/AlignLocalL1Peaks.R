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
BamFilePath     <- 'D:/L1polymORF/Data/BZ_NA12878L1capt5-9kb_subreads_hg19.bam'

# Specify the minimum read depth to create alignment
MinReadDepth <- 5

#######                       
# Get L1 ranges                    
#######                                     

cat("Getting reference L1 ranges \n")

# Read in table with known L1 
L1Catalogue <- read.csv(L1CataloguePath, as.is = T)

# Retain only entries with mapped L1 insertions and allele 1
blnL1Mapped       <- !is.na(L1Catalogue$start_HG38) 
blnAllele1        <- L1Catalogue$Allele == 1 
blnInRef          <- L1Catalogue$end_HG38 - L1Catalogue$start_HG38 > 5900
L1CatalogL1Mapped <- L1Catalogue[which(blnL1Mapped & blnAllele1 & blnInRef),]
sum(blnL1Mapped)
sum(blnAllele1, na.rm = T)
sum(blnInRef, na.rm = T)

# Lift coordinates and get genomic ranges for catalog L1 on hg19
LiftOverList <- LiftoverL1Catalog(L1CatalogL1Mapped,
    ChainFilePath = ChainFile)
L1CatalogGR <- LiftOverList$GRCatalogue_hg19
L1CatalogL1Mapped <- L1CatalogL1Mapped[LiftOverList$idxUniqueMapped, ]
width(L1CatalogGR)

#####                                   
# Analyze insertion flanks                       
#####

# Initialize a result table
ResultTable <- L1CatalogL1Mapped[,c("Accession", "Chromosome", 
        "strand_L1toRef", "start_HG38", "end_HG38")]
ResultTable$Coverage5P  <- NA
ResultTable$Coverage3P  <- NA
ResultTable$ReadWidth5P <- NA
ResultTable$ReadWidth3P <- NA

# Create a vector for start and end of a range
StartEnd <- c("start", "end")
Rshift    <- c(-10, 10)

# Loop through L1s and create local alignments
FastaFileNames <- c()
for (i in 1:nrow(ResultTable)){
  Acc <- ResultTable$Accession[i]
  Str <- ResultTable$strand_L1toRef[i]
  cat("Analyzing", Acc, "\n")
  for (j in 1:length(StartEnd)){
    if ((Str == "+" & j == 1) |  (Str == "-" & j == 2)){
      L1Side <- '5P'
    } else {
      L1Side <- '3P'
    }
    
    # Get range to intersect with reads (shifted to the out side of L1)
    R <- resize(L1CatalogGR[i], 1, fix = StartEnd[j], ignore.strand = T)
    strand(R) <- "*"
    R <- shift(R, Rshift[j])
    paramScan  <- ScanBamParam(which = R, what = scanBamWhat())
    Reads <- scanBam(BamFilePath, param = paramScan)
    
    # Retain only reads that actually intersect (bug in scanBam?)
    if (length(Reads[[1]]$pos) > 0){
      ReadGR <- GRanges(seqnames = Reads[[1]]$rname[1],
                        ranges = IRanges(Reads[[1]]$pos, width = width(Reads[[1]]$seq)))
      blnOverlap <- overlapsAny(ReadGR, R)
      Reads[[1]] <- lapply(Reads[[1]], function(x) x[blnOverlap])
    }
    
    # Update coverage at flank
    NrReadCol <- paste('Coverage', L1Side, sep = "")
    ResultTable[i, NrReadCol] <- length(Reads[[1]]$pos)
    
    cat("Start", start(R), "\n")
    cat("Strand", Str, "\n")
    cat("The coverage at", i, NrReadCol, "equals", length(Reads[[1]]$pos), "\n")
    if (length(Reads[[1]]$pos) > MinReadDepth){
      
      # Get common start and end for all reads
      StartAll <- max(Reads[[1]]$pos)
      EndAll   <- min(Reads[[1]]$pos + width(Reads[[1]]$seq))
      RelStart <- StartAll - Reads[[1]]$pos + 1
      RelEnd   <- RelStart + EndAll - StartAll - 1
      
      # Get sequences in common range 
      Seqs <- Reads[[1]]$seq
      cat("Getting sequences\n")
      for (k in 1:length(Seqs)){
        cat("Sequence", k, "with strand", Reads[[1]]$strand[k], "\n")
        #if (i == 20 & k == 7 & L1Side == '5P') browser()
        if (Reads[[1]]$strand[k] == 2){
          
          Seqs[k] <- reverseComplement(Seqs[k])
        }
        Seqs[k] <- as(Seqs[k][[1]][RelStart[k]:RelEnd[k]], "DNAStringSet")
      }
      
      # Get the reference sequence in common range
      RefSeq      <- getSeq(BSgenome.Hsapiens.UCSC.hg19, 
                            names = Reads[[1]]$rname[1], 
                       start = StartAll, end = EndAll)
      SeqsWithRef  <- c(as(RefSeq, "DNAStringSet"), Seqs)
      ReadWidthCol <- paste("ReadWidth", L1Side, sep = "")
      ResultTable[i, ReadWidthCol] <- width(SeqsWithRef)[1]
      
      # Create a file name for fasta file
      FastaFileName <- paste(FastaFolder, "FlankSeqs_", Acc, "_", L1Side, ".fas",
                             sep = "")
      cat("Writing sequences to", FastaFileName, "\n")
      writeFasta(SeqsWithRef, FastaFileName)
      FastaFileNames <- c(FastaFileNames, FastaFileName)
    }
  }
}

# Loop through fasta files and align them
for (FastaFile in FastaFileNames){
  AlignedFile <- gsub(".fas", "_aligned.fas", FastaFile)
  cat("Aligning", FastaFile, "\n")
  run_MUSCLE(InputPath = FastaFile, OutputPath = AlignedFile) 
}