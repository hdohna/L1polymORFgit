# The script below reads a file with sample info from the 1000 genome data and
# creates one file per population to be used with vcf tools

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Specify file paths
G1000SamplePath <- 'D:/L1polymORF/Data/1000GenomeSampleInfo.txt'
SuperPopPath    <- 'D:/L1polymORF/Data/1000G_SuperPop_'
PopPath         <- 'D:/L1polymORF/Data/1000G_Pop_'
CmdPathSuperPops <- 'D:/L1polymORFgit/Scripts/VcfFstCmds_SuperPops'
CmdPathPops   <- 'D:/L1polymORFgit/Scripts/VcfFstCmds_Pops'

# Read in table with Info about 1000 genome samples 
SampleInfo_1000Genome <- read.table(G1000SamplePath, as.is = T, header = T)

# Get names of super-populations
SuperPops <- unique(SampleInfo_1000Genome$super_pop)
SuperPopFiles <- c()
for (Pop in SuperPops){
  SuperPopFile <- paste(SuperPopPath, Pop, sep = "")
  Samples <- SampleInfo_1000Genome$sample[SampleInfo_1000Genome$super_pop == Pop]
  writeLines(Samples, SuperPopFile)
  SuperPopFiles <- c(SuperPopFiles, SuperPopFile)
}
substr(SuperPopFiles, 20, nchar(SuperPopFiles))
# Create vcftools commands for Fst
writeLines( c('vcftools --vcf L1all.vcf \\',
              paste('--weir-fst-pop', substr(SuperPopFiles, 20, nchar(SuperPopFiles)), '\\'),
              '--out L1all.weir.fst.SuperPops'),
            CmdPathSuperPops)

# Get names of super-populations
Pops <- unique(SampleInfo_1000Genome$pop)
PopFiles <- c()
for (Pop in Pops){
  PopFile <- paste(PopPath, Pop, sep = "")
  Samples <- SampleInfo_1000Genome$sample[SampleInfo_1000Genome$pop == Pop]
  writeLines(Samples, PopFile)
  PopFiles <- c(PopFiles, PopFile)
}

# Create vcftools commands for Fst
writeLines( c('vcftools --vcf L1all.vcf \\',
   paste('--weir-fst-pop', substr(PopFiles, 20, nchar(PopFiles)), '\\'),
   '--out L1all.weir.fst.Pops'),
   CmdPathPops)

