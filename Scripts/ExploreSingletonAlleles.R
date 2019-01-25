# The following script explores why the singleton alleles are not always 1 and 3 

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load 1000 genome data
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

for (Chr in c(6:10)) {
  cat("********   Analyzing chromosome", Chr, "    **********\n")
  
  # Read in singletons for current chromosome
  SingletonPath  <- paste("D:/L1polymORF/Data/Singleton_SNP_chr", Chr, sep = "")
  Singletons     <- read.table(SingletonPath)
  SCols          <- GetSingletonColumns(Singletons)
  
  cat("Unique alleles:", unique(SCols$Allele), "\n")
#  plot(SCols$Allele)
}


# Explore chromosome 1
Singletons <- read.table("D:/L1polymORF/Data/Singleton_SNP_chr1")
SCols      <- GetSingletonColumns(Singletons)
idxMinWrongAllele <- min(which(SCols$Allele %in% c(0, 2)))
Singletons$V2[idxMinWrongAllele + 1] 

# Read vcf file
Line1Vcf  <- read.table("D:/L1polymORF/Data/LINE1chr1.vcf", as.is = T)

L1_1000G[L1_1000G$POS == L1_1000G, SampleColumns[SCols$Col[idxMinWrongAllele]]

# Read in singletons created by vcf tools
Chr1Singl <- read.table("D:/L1polymORF/Data/chr1.singletons", header = T)
