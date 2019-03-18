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

#############################################
#                                           #
#   Analyze L1 detection from standard      #
#          simulated files                  #
#                                           #
#############################################

# Specify simulation directory
SimDir <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT_single/"

# Get names of vcf files
VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/", full.names = F)
VcfDirs <- VcfDirs[VcfDirs %in% paste(SampleColumns, "fullGenome", sep = "_")]
VcfDirs[1:10]
VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/", VcfDirs, 
                  "/LINE1.final_comp.vcf", sep = "")
FileSizes <- file.info(VcfFiles)$size
VcfFiles[1:10]
sum(is.na(FileSizes))
sum(FileSizes == 0, na.rm = T)
sum(!file.exists(VcfFiles))

# Read in vcf file and create genomic ranges
L1Vcf <- NULL
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
