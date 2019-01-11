# The script below identifies deletions that are likely L1 polymorphis,s

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Specify file paths
DataFolder            <- 'D:/L1polymORF/Data/'
MeltInsPath         <- "D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf"
MeltDelPath         <- "D:/L1polymORF/Data/DEL.final_comp.vcf"
ChrLPath            <- 'D:/L1polymORF/Data/ChromLengthsHg19.Rdata'
InputPath           <- 'D:/L1polymORF/Data/SingletonAnalysis_unphased.RData'
L1RefPath           <- 'D:/L1polymORF/Data/L1HS_repeat_table_Hg19.csv'
L1RefRangePath      <- 'D:/L1polymORF/Data/L1RefRanges_hg19.Rdata'
RegrOutputPath      <- "D:/L1polymORF/Data/L1RegressionResults.RData"
SelectTabOutPath    <- "D:/L1polymORF/Data/L1SelectionResults_MELT.csv"
SelectGenTabOutPath <- "D:/L1polymORF/Data/L1SelectionGeneResults_MELT.csv"
SelectResultOutPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT.RData"
SelectWithinGenTabOutPath <- "D:/L1polymORF/Data/L1SelectionWithinGeneResults_MELT.csv"
SelectSingletonTabOutPath <- "D:/L1polymORF/Data/L1SelectionSingletonResults_MELT.csv"

# Read in ranges of non-polymorphic
L1NoPolyGR <- import.bed("D:/L1polymORF/Data/L1GRangesNoPoly.bed")


# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(DataFolder, pattern = "_L1NoPoly_Indels", 
                       full.names = T)
InFile <- AllFiles[1]
StartV <- NULL
EndV   <- NULL
ChromV <- NULL
for (InFile in AllFiles){
  VcfFile <- read.table(InFile, stringsAsFactors = F)
  Chrom   <- strsplit(InFile, "_")[[1]][2]
  cat("Analyzing", Chrom, ":\n")
  idxCNV  <- which(VcfFile$V5 == "<CN0>")
  StartV  <- c(StartV, VcfFile$V2[idxCNV])
  ChromV  <- c(ChromV, rep(Chrom, length(idxCNV)))
  Ends <- sapply(VcfFile$V8[idxCNV], function(x){
    Split1 <- strsplit(x, ";")[[1]]
    EndPart <- grep("END=", Split1, value = T)
    if(length(grep("CIEND=", EndPart)) > 0){
      EndPart <- EndPart[-grep("CIEND=", EndPart)]
    }
    if (length(EndPart) > 0){
      as.numeric(strsplit(EndPart, "END=")[[1]][2])
    } else {
      NA
    }
    
  })
  EndV <- c(EndV, unlist(Ends))
  cat("Added", length(idxCNV), "start and chromosome values and", length(Ends),
      "end values\n\n")
  
}
EndV - StartV

# Create genomic ranges for start and end of variants
StartGR <- GRanges(seqnames = substr(ChromV, 4, nchar(ChromV)), IRanges(start = StartV - 10, end = StartV + 10))
EndGR   <- GRanges(seqnames = substr(ChromV, 4, nchar(ChromV)), IRanges(start = EndV - 10, end = EndV + 10))

# Get indicator for L1s that overlap with either end of the CNV
blnOL   <- overlapsAny(L1NoPolyGR, StartGR) | overlapsAny(L1NoPolyGR, EndGR)

# Check whether probability of overlap depends on lend and full-length indicator
blnFull <- width(L1NoPolyGR) >= 6000
boxplot(width(L1NoPolyGR) ~ blnOL)
GLM_OL <- glm(blnOL ~ width(L1NoPolyGR) + blnFull, family = binomial)
summary(GLM)
