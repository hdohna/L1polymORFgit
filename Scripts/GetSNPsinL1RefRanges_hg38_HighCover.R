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
source('/home/hb54/L1polymORFgit/Scripts/_Start_L1polymORF_AUB.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

ChrLPath        <- '/home/hb54/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName <- "/home/hb54/RefSeqData/L1HS_repeat_table_Hg19.csv"
BedPath_L1      <- '/home/hb54/RefSeqData/L1HSRefRanges_hg38.bed'
BedPath_L1_left      <- '/home/hb54/RefSeqData/L1HSRefRanges_hg38_left1000.bed'
BedPath_L1_right      <- '/home/hb54/RefSeqData/L1HSRefRanges_hg38_right1000.bed'
Folder1000GenomesIn  <- "/scratch/hb54/1000G_high_coverage/"
Folder1000GenomesOut <- '/home/hb54/1000GenomeData/'
FullLength <- 6000

#######################################
#                                     #
#      Get SNPs per chromosome         #
#                                     #
#######################################

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(Folder1000GenomesIn, pattern = "CCDG_13607_B01_GRM_WGS_2019-02-19", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.gz.tbi", AllFiles)]
AllFiles <- AllFiles[-grep("others.recalibrated", AllFiles)]
AllFiles
CurrentFile <- AllFiles[1]

# Initialize vector of output files 
OutFilesL1 <- c()
OutFilesL1left <- c()
OutFilesL1right <- c()

# Loop over all chromosomes and get snps on that chromosome
for (CurrentFile in  AllFiles){
  
  Chrom <- strsplit(CurrentFile, "19_")[[1]][2]
  Chrom <- strsplit(Chrom, "\\.")[[1]][1]
  
  # SNPs in L1s
  ScriptFile <- paste0("runVcftools_L1_", Chrom)
  ScriptName <- paste0("Vcftools_L1_", Chrom)
  OutPath <- paste0(Folder1000GenomesOut, "SNPsInL1_HighCover_", Chrom)
  OutFilesL1 <- c(OutFilesL1, OutPath)
  
  # Run vcftools script for each set of ranges 
  CreateAndCallSlurmScript(file = ScriptFile, 
                           scriptName = ScriptName,
                           SlurmCommandLines = c("module load vcftools",
                                                 paste("vcftools --gzvcf", CurrentFile, 
                                                       "--bed", BedPath_L1, 
                                                       "--recode --recode-INFO-all",
                                                       "--out", OutPath)
                           ),
                           RunTime = '12:00:00',
                           Mem = '50G')

  # SNPs in left flank of L1s
  ScriptFile <- paste0("runVcftools_L1_left_", Chrom)
  ScriptName <- paste0("Vcftools_L1_left_", Chrom)
  OutPath <- paste0(Folder1000GenomesOut, "SNPsInL1left_HighCover_", Chrom)
  OutFilesL1left <- c(OutFilesL1left, OutPath)
  
  # Run vcftools script for each set of ranges 
  CreateAndCallSlurmScript(file = ScriptFile, 
                           scriptName = ScriptName,
                           SlurmCommandLines = c("module load vcftools",
                                                 paste("vcftools --gzvcf", CurrentFile, 
                                                       "--bed", BedPath_L1_left, 
                                                       "--recode --recode-INFO-all",
                                                       "--out", OutPath)
                           ),
                           RunTime = '12:00:00',
                           Mem = '50G')
  
  # SNPs in right flank of L1s
  ScriptFile <- paste0("runVcftools_L1_right_", Chrom)
  ScriptName <- paste0("Vcftools_L1_right_", Chrom)
  OutPath <- paste0(Folder1000GenomesOut, "SNPsInL1right_HighCover_", Chrom)
  OutFilesL1right <- c(OutFilesL1right, OutPath)
  
  # Run vcftools script for each set of ranges 
  CreateAndCallSlurmScript(file = ScriptFile, 
                           scriptName = ScriptName,
                           SlurmCommandLines = c("module load vcftools",
                                                 paste("vcftools --gzvcf", CurrentFile, 
                                                       "--bed", BedPath_L1_right, 
                                                       "--recode --recode-INFO-all",
                                                       "--out", OutPath)
                           ),
                           RunTime = '12:00:00',
                           Mem = '50G')
  
}

# Concatenate different vcf files for SNPs in L1 
ConcatVCF <-  data.frame()
OutFilesVCf_L1 <- list.files(Folder1000GenomesOut, pattern = "SNPsInL1_HighCover_", 
                             full.names = T)
for (File in OutFilesVCf_L1){
  NewFile <- ReadVCF(File)
  ConcatVCF <- rbind(ConcatVCF, NewFile)
}
write.table(ConcatVCF, file = "/home/hb54/1000GenomeData/SNPsInL1_highCover_all.vcf")

# Concatenate different vcf files for SNPs in L1 left flank 
ConcatVCF <-  data.frame()
OutFilesVCf_L1left <- list.files(Folder1000GenomesOut, pattern = "SNPsInL1left_HighCover_", 
                             full.names = T)
OutFilesVCf_L1left
for (File in OutFilesVCf_L1left){
  NewFile <- ReadVCF(File)
  ConcatVCF <- rbind(ConcatVCF, NewFile)
}
write.table(ConcatVCF, file = "/home/hb54/1000GenomeData/SNPsInL1left_highCover_all.vcf")

# Concatenate different vcf files for SNPs in L1 right flank 
ConcatVCF <-  data.frame()
OutFilesVCf_L1right <- list.files(Folder1000GenomesOut, pattern = "SNPsInL1right_HighCover_", 
                                 full.names = T)
OutFilesVCf_L1right
for (File in OutFilesVCf_L1right){
  NewFile <- ReadVCF(File)
  ConcatVCF <- rbind(ConcatVCF, NewFile)
}
write.table(ConcatVCF, file = "/home/hb54/1000GenomeData/SNPsInL1right_highCover_all.vcf")
