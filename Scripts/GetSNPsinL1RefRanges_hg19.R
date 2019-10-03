##############################################
#
# General description:
#
#   The following script reads the repeat masker table for hg19
#   and creates genomic ranges of L1HS, creates bed files and
#   gets the number of SNPs in these ranges using Vcftools

# Input:
#

# Output:
#   
#    : ...

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

ChrLPath           <- '/labs/dflev/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName    <- "/labs/dflev/hzudohna/RefSeqData/L1HS_repeat_table_Hg19.csv"
OutBedPath_L1      <- '/labs/dflev/hzudohna/RefSeqData/L1HSRefRanges_hg19.bed'
OutBedPath_L1Left  <- '/labs/dflev/hzudohna/RefSeqData/L1HSRef_Left1000_hg19.bed'
OutBedPath_L1Right <- '/labs/dflev/hzudohna/RefSeqData/L1HSRef_Right1000_hg19.bed'
OutVcfPath_L1      <- '/labs/dflev/hzudohna/1000Genomes/VariantsInL1'
OutVcfPath_L1Left  <- '/labs/dflev/hzudohna/1000Genomes/VariantsInL1_leftFlank'
OutVcfPath_L1Right <- '/labs/dflev/hzudohna/1000Genomes/VariantsInL1_rightFlank'

FullLength <- 6000
FlankSize  <- 1000

#######################################
#                                     #
#    Read in data                     #
#                                     #
#######################################

# Load vector with chromosome lengths
load(ChrLPath)

# Read in table with L1 ranges
L1Table <- read.csv(L1TableFileName, as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))
L1Table$genoStartMinus1000 <- L1Table$genoStart - FlankSize
L1Table$genoEndPlus1000    <- L1Table$genoEnd + FlankSize

# Create names of subfamilies
repNChar    <- as.character(L1Table$repName)
SubFamiliesLookUp <- sapply(unique(repNChar), function(x){
  c(Name = x,
    SubFam = paste(strsplit(x, '[1-9]')[[1]][1:2], collapse = "1"))
})
NameMatch   <- match(repNChar, SubFamiliesLookUp["Name",])
SubFamilies <- SubFamiliesLookUp["SubFam", NameMatch]

# Create GRanges objects with L1 Seqences
L1GRanges <- makeGRangesFromDataFrame(L1Table, seqnames.field = "ChrNr",
                                      start.field = "genoStart",
                                      end.field = "genoEnd")

# Get flanking ranges
L1LeftFlankRanges <- makeGRangesFromDataFrame(L1Table, 
                                              seqnames.field = "ChrNr",
                                              start.field = "genoStartMinus1000",
                                              end.field = "genoStart")

L1RightFlankRanges <- makeGRangesFromDataFrame(L1Table, 
                                              seqnames.field = "ChrNr",
                                              start.field = "genoEnd",
                                              end.field = "genoEndPlus1000")


#######################################
#                                     #
#    Export data and run vcftools     #
#                                     #
#######################################

# Export bed ranges
export.bed(L1GRanges, con = OutBedPath_L1)
export.bed(L1LeftFlankRanges, con = OutBedPath_L1Left)
export.bed(L1RightFlankRanges, con = OutBedPath_L1Right)

# Run vcftools script for each set of ranges 
CreateAndCallSlurmScript(file = "runVcftools_L1", 
                         scriptName = "Vcftools_L1",
                         SlurmCommandLines = c("module load vcftools",
                            paste("vcftools --gzvcf  /reference/1KGenomes/Variants/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz --bed",
                                  OutBedPath_L1, "--recode --remove-indels --out",
                                  OutVcfPath_L1)
                         ),
                         RunTime = '12:00:00',
                         Mem = '200G')

CreateAndCallSlurmScript(file = "runVcftools_L1left", 
                         scriptName = "Vcftools_L1left",
                         SlurmCommandLines = 
                           c("module load vcftools",
                              paste("vcftools --gzvcf  /reference/1KGenomes/Variants/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz --bed",
                                    OutBedPath_L1Left, "--recode --remove-indels --out",
                                    OutVcfPath_L1Left)),
                         RunTime = '12:00:00',
                         Mem = '200G')

CreateAndCallSlurmScript(file = "runVcftools_L1Right", 
                         scriptName = "Vcftools_L1Right",
                         SlurmCommandLines = 
                           c("module load vcftools",
                             paste("vcftools --gzvcf  /reference/1KGenomes/Variants/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz --bed",
                                   OutBedPath_L1Right, "--recode --remove-indels --out",
                                   OutVcfPath_L1Right)),
                         RunTime = '12:00:00',
                         Mem = '200G')

#######################################
#                                     #
#      Calculate HWE p-values         #
#                                     #
#######################################

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
DataFolder <- "/labs/dflev/hzudohna/1000Genomes/"
AllFiles <- list.files(DataFolder, pattern = "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.", AllFiles)]
PopFiles <- list.files(DataFolder, pattern = "1000G_SuperPop", 
                       full.names = T)
HWEOutFilesL1    <- c()
HWEOutFilesLeft  <- c()
HWEOutFilesRight <- c()
for (InFile in AllFiles){
  InFileSplit  <- strsplit(InFile, "\\.")[[1]]
  OutFile      <- paste(InFileSplit[1:2], collapse = "")
  OutFileLeft  <- gsub("ALL", "HWE_LeftFlank_", OutFile)
  OutFileRight <- gsub("ALL", "HWE_RightFlank_", OutFile)
  cat("\n*****   Calculating HWE for", InFileSplit[2], "    ****************\n\n")
  for (PopFile in PopFiles){
    PopFileSplit  <- strsplit(PopFile, "\\_")[[1]]
    OutFileL1     <- gsub("ALL", paste("HWE_L1_Pop", PopFileSplit[3], sep = "_"), OutFile)
    cat("HWE for", PopFileSplit[3], "\n")
    FileName <- paste("vcfHWE", InFileSplit[2], PopFileSplit[3], sep = "_")
    CreateAndCallSlurmScript(file = FileName, 
                             scriptName = FileName,
                             SlurmCommandLines = 
                               c("module load vcftools",
                                 paste("vcftools --vcf", InFile, 
                                       "--bed", OutBedPath_L1, 
                                       "--keep", PopFile,
                                       "--hardy --out",
                                       OutFileL1)),
                             RunTime = '12:00:00',
                             Mem = '200G')
    
  }
}

# Get all files with HWE for L1 and concatenate them
HWEOutFilesL1 <- list.files(DataFolder, pattern = "HWE_L1_Pop", 
                            full.names = T)
# Make genomic ranges
HWEGR <- GRanges()

for (i in 1:length(HWEOutFilesL1)){
  
  cat("Processing file", i, "of", length(HWEOutFilesL1), "\n")
  # Read in HWE data
  HWE_Local  <- read.table(HWEOutFilesL1[i], header = T)
  
  # Get observed and expected heterozygotes
  blnObsHetGreater <- as.data.frame(t(sapply(1:nrow(HWE_Local), function(i){
    ObsHetChar <- strsplit(as.character(HWE_Local$OBS.HOM1.HET.HOM2.[i]), "/")[[1]][2]
    ExpHetChar <- strsplit(as.character(HWE_Local$E.HOM1.HET.HOM2.[i]), "/")[[1]][2]
    as.numeric(ObsHetChar) > as.numeric(ExpHetChar)
  })))
  
  # Subset to get only entries with P <= 0.05 and observed heterozygosity
  # greater 
  HWE_Local <- HWE_Local[which(HWE_Local$P <= 0.05 & blnObsHetGreater), ]
  
  if(nrow(HWE_Local) > 0){
    # Get genomic ranges
    HWE_Local$chrom <- paste("chr", HWE_Local$CHR, sep = "")
    HWEGR_Local <- makeGRangesFromDataFrame(HWE_Local, seqnames.field = "chrom",
                                            start.field = "POS",
                                            end.field = "POS")
    HWEGR <- union(HWEGR, HWEGR_Local)
  }
  
  cat("Total SNPs found:", length(HWEGR), "\n")
}


Explanation <- c("This R workspace contains genomic ranges of SNPs that",
                 "differ significantly from HWE and have more heterozygotes",
                 "than expected.")
save(list = c("HWEGR", "Explanation"), 
     file = paste(DataFolder, "SNPsDifferFromHWE.RData", sep = ""))

# # Get all files with HWE for left flank and concatenate them
# HWEOutFilesLeft <- list.files(DataFolder, pattern = "HWE_Left", 
#                             full.names = T)
# HWE_Left <- read.table(HWEOutFilesLeft[1], header = T)
# for (i in 2:length(HWEOutFilesLeft)){
#   HWE_Local <- read.table(HWEOutFilesLeft[i], header = T)
#   HWE_Left <- rbind(HWE_Left, HWE_Local)
# }
# 
# # Get all files with HWE for rigth flank and concatenate them
# HWEOutFilesRight <- list.files(DataFolder, pattern = "HWE_Right", 
#                               full.names = T)
# HWE_Right <- read.table(HWEOutFilesRight[1], header = T)
# for (i in 2:length(HWEOutFilesRight)){
#   HWE_Local <- read.table(HWEOutFilesRight[i], header = T)
#   HWE_Right <- rbind(HWE_Right, HWE_Local)
# }

Explanation <- c("This R workspace contains the following dataframes",
                 "HWE_L1: Hardy Weinberg equilibrium of SNPs within reference L1",
                 "HWE_Left: Hardy Weinberg equilibrium of SNPs 1000 bp left of reference L1",
                 "HWE_Right: Hardy Weinberg equilibrium of SNPs 1000 bp right of reference L1")
save(list = c("HWE_L1", "HWE_Left", "HWE_Right", "Explanation"), 
     file = paste(DataFolder, "HWE_L1_and Flanks.RData", sep = ""))
                 
