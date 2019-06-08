# The following script runs commands pair read positions on L1 and genome.

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load necessary packages
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)

# Path to reference data and for 1000 genomes
PathL1simulated <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"
OutSuffix       <- "LINE1.hum.readpositions"

# Load necessary objects
load('/labs/dflev/hzudohna/1000Genomes/GRanges_L1_1000Genomes.RData')

# Get all files with simulated genomes
SimDirs <- paste(PathL1simulated, SampleColumns, sep = "")
SimDirs <- SimDirs[dir.exists(SimDirs)]
SimDirs

##################################################################
#                                                                #                             
#         Construct and run commands to pair read positions      #
#                                                                #                             
##################################################################

# Loop over directories 
OutFiles <- NULL
for (SimDir in SimDirs){
  DirSplit <- strsplit(SimDir, "/")[[1]]
  SampleID <- DirSplit[length(DirSplit)]
  Line1Bam  <- paste(SimDir, "/hg19_", SampleID, 
                     "_sorted.LINE1.aligned.final.sorted.bam", sep = "")
  GenomeBam  <- paste(SimDir, "/LINE1.merged.hum_breaks.sorted.bam", sep = "")
  OutputFile <- paste(SimDir, "/", OutSuffix, sep = "")
  Cmd <- paste("sbatch sb_PairL1GenomeReads", Line1Bam, GenomeBam, OutputFile)
  system(Cmd)
  OutFiles <- c(OutFiles, OutputFile)
}

##################################################################
#                                                                #                             
#       Match read positions to insertions                       #
#                                                                #                             
##################################################################

VcfFiles <- paste(SimDirs, "/LINE1.final_comp.vcf", sep = "")

# Loop over directories
SimDir <- SimDirs[2]
L1Detect <- NULL
for (SimDir in SimDirs[file.exists(VcfFiles)]){

  DirSplit <- strsplit(SimDir, "/")[[1]]
  SampleID <- DirSplit[length(DirSplit)]
  cat("Processing", SampleID, "\n")
  VcfFile    <- paste(SimDir, "LINE1.final_comp.vcf", sep = "/")
  OutputFile <- paste(SimDir, "/", OutSuffix, sep = "")
  
  # Read in vcf file and make genomic ranges
  VcfFile <- ReadVCF(VcfFile)
  VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                      start.field = "POS",
                                      end.field = "POS")
  VcfGR <- resize(VcfGR, width = 100, fix = "center")

  # Read in position pair file and make genomic ranges
  PosPairFile <- read.table(OutputFile,header = T)
  PosPairGR   <- makeGRangesFromDataFrame(PosPairFile, seqnames.field = "CHROM",
                                      start.field = "POS",
                                      end.field = "POS")

  # Create overlaps and get minimum and maximum on L1 per overlap
  OLSimPosPair <- findOverlaps(VcfGR, PosPairGR) 
  MinL1Pos <- aggregate(PosPairFile$POSL1[OLSimPosPair@to], 
                        by = list(OLSimPosPair@from),
                        FUN = min)
  colnames(MinL1Pos) <- c("idxVcf", "MinL1")
  MaxL1Pos <- aggregate(PosPairFile$POSL1[OLSimPosPair@to], 
                        by = list(OLSimPosPair@from),
                        FUN = max)
  colnames(MaxL1Pos) <- c("idxVcf", "MaxL1")
  L1MinMax <- merge(MinL1Pos, MaxL1Pos)
  head(L1MinMax)
  
  # Get estimated L1 length and genotype from vcf file
  L1widthVcf <- sapply(VcfFile$INFO, GetFromVcfINFO_SVLength)
  L1GenoVcf  <- sapply(VcfFile[,ncol(VcfFile)], GetFromVcfGeno_GenoNum)
  L1StartEndVcf  <- t(sapply(VcfFile$INFO, GetFromVcfINFO_MELT_L1StartEnd))
  L1StartEndVcf
  
  # Get sample ID from vcf column and get index of L1 insertions that occur in
  # that sample
  SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
  idxL1    <- which(L1_1000G[,SampleID] > 0)
  SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
                                       start.field = "POS",
                                       end.field = "POS")
  
  # Determine overlaps between simulated and detected L1
  OLSimDetect <- findOverlaps(VcfGR, SampleGR)
  PosMatch    <- match(OLSimDetect@from, L1MinMax$idxVcf)
  
  # Create data.frame that keeps track of L1 present in 1000 Genomes and 
  # their detection status
  L1DetectNew <- data.frame(Chrom  = L1_1000G$CHROM[idxL1],
                            PosTrue  = L1_1000G$POS[idxL1],
                            PosEst      = NA,
                            blnDetect   = overlapsAny(SampleGR, VcfGR),
                            L1GenoTrue  = L1_1000G[idxL1, SampleID],
                            L1widthTrue = L1_1000G$InsLength[idxL1],
                            L1StartTrue = L1_1000G$L1Start[idxL1],
                            L1EndTrue   = L1_1000G$L1End[idxL1],
                            L1widthEst  = NA,
                            L1StartEst  = NA,
                            L1EndEst  = NA,
                            L1StartEst_DiscReads  = NA,
                            L1EndEst_DiscReads  = NA,
                            L1GenoEst   = NA,
                            EstFilter   = NA,
                            SampleID = SampleID)
  L1DetectNew$PosEst[OLSimDetect@to]     <- VcfFile$POS[OLSimDetect@from]
  L1DetectNew$L1widthEst[OLSimDetect@to] <- L1widthVcf[OLSimDetect@from]
  L1DetectNew$L1StartEst[OLSimDetect@to] <- L1StartEndVcf[OLSimDetect@from, 1]
  L1DetectNew$L1widthEst[OLSimDetect@to] <- L1widthVcf[OLSimDetect@from]
  L1DetectNew$L1StartEst_DiscReads[OLSimDetect@to] <- L1MinMax$MinL1[PosMatch]
  L1DetectNew$L1EndEst_DiscReads[OLSimDetect@to]   <-  L1MinMax$MaxL1[PosMatch]
  L1DetectNew$L1GenoEst[OLSimDetect@to]  <- L1GenoVcf[OLSimDetect@from]
  L1DetectNew$EstFilter[OLSimDetect@to]  <- VcfFile$FILTER[OLSimDetect@from]
  L1Detect   <- rbind(L1Detect, L1DetectNew)
  
  # Create data.frame that keeps track of L1 not present in 1000 Genomes and 
  # their detection status
  blnPresent  <- overlapsAny(VcfGR, SampleGR)
  if(any(!blnPresent)){
    L1DetectNew <- data.frame(Chrom  = VcfFile$X.CHROM[!blnPresent],
                              PosTrue  = NA,
                              PosEst      = VcfFile$POS[!blnPresent],
                              blnDetect   = NA,
                              L1GenoTrue  = 0,
                              L1widthTrue = 0,
                              L1StartTrue = NA,
                              L1EndTrue   = NA,
                              L1widthEst  = L1widthVcf[!blnPresent],
                              L1StartEst  =  L1StartEndVcf[!blnPresent, 1],
                              L1StartEst_DiscReads = NA,
                              L1EndEst_DiscReads   = NA,
                              L1EndEst    =  L1StartEndVcf[!blnPresent, 2],
                              L1GenoEst   = L1GenoVcf[!blnPresent],
                              EstFilter   = VcfFile$FILTER[!blnPresent],
                              SampleID = SampleID)
    L1Detect   <- rbind(L1Detect, L1DetectNew)
  }
  
}
head(L1Detect)
