# The following script compares two different files that contain singletons. 
# One file is created by vcftools 

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load 1000 genome data
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

# Read in singletons created by vcf tools
Chr1SinglVcf <- read.table("D:/L1polymORF/Data/chr1.singletons", header = T)
Diff <- Chr1SinglVcf$POS[-1] - Chr1SinglVcf$POS[-nrow(Chr1SinglVcf)]
all(Diff >= 0)

# Read in singletons created by created by Singleton_awk_cmds
Singletons <- read.table("D:/L1polymORF/Data/Singleton_SNP_chr1")
SCols      <- GetSingletonColumns(Singletons)

# Match positions between files
PosMatch  <- match(Singletons$V2, Chr1SinglVcf$POS)
Chr1SinglVcf_matched <- Chr1SinglVcf[PosMatch,]

# Compare sample columns of file
all(SCols$Col == match(Chr1SinglVcf$INDV[PosMatch], SampleColumns))
blnNoMatch <- SCols$Col != match(Chr1SinglVcf$INDV[PosMatch], SampleColumns)
sum(blnNoMatch)
nrow(SCols)
table(SCols$Col[blnNoMatch])
