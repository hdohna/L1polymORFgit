# The following script gets vcf files of identified L1 insertions from
# 1000 genome by running MELT single option

# Source start script
#source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Load packages
library(GenomicRanges)

# Load 1000 genome data
#load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
load('/labs/dflev/hzudohna/1000Genomes/GRanges_L1_1000Genomes.RData')
SampleInfo <- read.table("/labs/dflev/hzudohna/1000Genomes/1000GenomeSampleInfo.txt", header = T)

#############################################
#                                           #
#   Concatenate vcf files                   #
#                                           #
#############################################

# Get names of vcf files
VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/", full.names = F)
VcfDirs <- VcfDirs[VcfDirs %in% paste(SampleColumns, "fullGenome", sep = "_")]
VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/", VcfDirs, 
                  "/LINE1.final_comp.vcf", sep = "")
FileSizes <- file.info(VcfFiles)$size
cat(sum(FileSizes == 0, na.rm = T), "empty vcf files\n")
cat(sum(!file.exists(VcfFiles)), "vcf files don't exist\n")

# Read in vcf file and create genomic ranges
L1Vcf <- NULL
VcfFile <- VcfFiles[1]
cat(sum(file.exists(VcfFiles)), "Vcf files exist\n")
cat(sum(FileSizes > 0), "Vcf files with data\n")
for (VcfFile in VcfFiles[which(file.exists(VcfFiles) & FileSizes > 0)]){
  
  ID <- strsplit(VcfFile, "/")[[1]][6]
  cat("Processing", ID, "\n")
  
  VcfFile <- ReadVCF(VcfFile)
  colnames(VcfFile)[10] <- "Genotype"
  VcfFile$SampleID <- ID
  L1Vcf   <- rbind(L1Vcf, VcfFile)
  
}

# Write out vcf file
write.table(L1Vcf, "/labs/dflev/hzudohna/1000Genomes/L1_SingleMELT_fullGenome_CombinedVcfs")

#############################################
#                                           #
#   Concatenate vcf files per               #
#       super-population                    #
#                                           #
#############################################

# Get all super populations
SuperPops <- unique(SampleInfo$super_pop)
unique(SampleInfo$pop)

for(SuperPop in SuperPops){
  
  cat("**********  Processing", SuperPop, "************* \n")
  
  # Get all sample columns of current population
  SampleColumns <- SampleInfo$sample[SampleInfo$super_pop = SuperPop]
  
  # Get names of vcf files
  VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/", full.names = F)
  VcfDirs <- VcfDirs[VcfDirs %in% paste(SampleColumns, "fullGenome", sep = "_")]
  VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/", VcfDirs, 
                    "/LINE1.final_comp.vcf", sep = "")
  FileSizes <- file.info(VcfFiles)$size
  cat(sum(FileSizes == 0, na.rm = T), "empty vcf files\n")
  cat(sum(!file.exists(VcfFiles)), "vcf files don't exist\n")
  
  # Read in vcf file and create genomic ranges
  L1Vcf <- NULL
  VcfFile <- VcfFiles[1]
  for (VcfFile in VcfFiles[which(file.exists(VcfFiles) & FileSizes > 0)]){
    
    ID <- strsplit(VcfFile, "/")[[1]][6]
    cat("Processing", ID, "\n")
    
    VcfFile <- ReadVCF(VcfFile)
    colnames(VcfFile)[10] <- "Genotype"
    VcfFile$SampleID <- ID
    L1Vcf   <- rbind(L1Vcf, VcfFile)
    
  }
  
  # Write out vcf file
  NewPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_SingleMELT_fullGenome_CombinedVcfs", SuperPop,
                   sep = "_")
  write.table(L1Vcf, NewPath)
  
}

#############################################
#                                           #
#   Concatenate vcf files per               #
#       population                    #
#                                           #
#############################################

# Get all super populations
SuperPops <- unique(SampleInfo$super_pop)
unique(SampleInfo$pop)

for(SuperPop in SuperPops){
  
  cat("**********  Processing", SuperPop, "************* \n")
  
  # Get all sample columns of current population
  SampleColumns <- SampleInfo$sample[SampleInfo$super_pop = SuperPop]
  
  # Get names of vcf files
  VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/", full.names = F)
  VcfDirs <- VcfDirs[VcfDirs %in% paste(SampleColumns, "fullGenome", sep = "_")]
  VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/", VcfDirs, 
                    "/LINE1.final_comp.vcf", sep = "")
  FileSizes <- file.info(VcfFiles)$size
  cat(sum(FileSizes == 0, na.rm = T), "empty vcf files\n")
  cat(sum(!file.exists(VcfFiles)), "vcf files don't exist\n")
  
  # Read in vcf file and create genomic ranges
  L1Vcf <- NULL
  VcfFile <- VcfFiles[1]
  for (VcfFile in VcfFiles[which(file.exists(VcfFiles) & FileSizes > 0)]){
    
    ID <- strsplit(VcfFile, "/")[[1]][6]
    cat("Processing", ID, "\n")
    
    VcfFile <- ReadVCF(VcfFile)
    colnames(VcfFile)[10] <- "Genotype"
    VcfFile$SampleID <- ID
    L1Vcf   <- rbind(L1Vcf, VcfFile)
    
  }
  
  # Write out vcf file
  NewPath <- paste("/labs/dflev/hzudohna/1000Genomes/L1_SingleMELT_fullGenome_CombinedVcfs", SuperPop,
                   sep = "_")
  write.table(L1Vcf, NewPath)
  
}

