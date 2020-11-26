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
#library(rtracklayer)

ChrLPath        <- '/home/hb54/RefSeqData/ChromLengthsHg19.Rdata'
L1TableFileName <- "/home/hb54/RefSeqData/L1HS_repeat_table_Hg19.csv"
BedPath_L1      <- '/home/hb54/RefSeqData/L1HSRefRanges_hg19.bed'
Folder1000GenomesIn  <- "/scratch/shared/projects/biol_290az_390af/1000GenomeData/"
Folder1000GenomesOut <- '/home/hb54/1000GenomeData/'
FullLength <- 6000

#######################################
#                                     #
#      Get SNPs per chromosome         #
#                                     #
#######################################

# Get file names, loop over files and do the filtering
# Example file name: ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
AllFiles <- list.files(Folder1000GenomesIn, pattern = "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 
                       full.names = T)
AllFiles <- AllFiles[-grep("vcf.gz.tbi", AllFiles)]
AllFiles
CurrentFile <- AllFiles[1]

# Initialize vector of output files 
OutFiles <- c()

# Loop over all chromosomes and get snps on that chromosome
for (CurrentFile in  AllFiles){
  
  Chrom <- strsplit(CurrentFile, "\\.")[[1]][2]
  ScriptFile <- paste0("runVcftools_L1_", Chrom)
  ScriptName <- paste0("Vcftools_L1_", Chrom)
  OutPath <- paste0(Folder1000GenomesOut, "SNPsInHG002_", Chrom)
  OutFiles <- c(OutFiles, OutPath)
  
  # Run vcftools script for each set of ranges 
  CreateAndCallSlurmScript(file = "runVcftools_L1", 
                           scriptName = "Vcftools_L1",
                           SlurmCommandLines = c("module load vcftools",
                                                 paste("vcftools --gzvcf", CurrentFile, 
                                                       "--bed", BedPath_L1, 
                                                       "--recode --recode-INFO-all --indv HG002",
                                                       "--out", OutPath)
                           ),
                           RunTime = '12:00:00',
                           Mem = '200G')
  
}

# Concatenate all SNP files
# OutFilesVCf <- paste0(OutFiles, ".recode.vcf")
# ConcatCmd <- paste("cat", paste(OutFilesVCf, collapse = " "), "> /home/hb54/1000GenomeData/SNPsInHG002_all")
# system(ConcatCmd)

# Concatenate different vcf files 
ConcatVCF <-  data.frame()
for (File in OutFilesVCf){
  NewFile <- ReadVCF(File)
  ConcatVCF <- rbind(ConcatVCF, NewFile)
}
write.table(ConcatVCF, file = "/home/hb54/1000GenomeData/SNPsInHG002_all.vcf")
