library("VariantAnnotation")
library("VariantTools")
library("seqinr")

source("https://bioconductor.org/biocLite.R")
biocLite("VariantTools")

GenomeList <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
Genome <- paste(GenomeList[[1]], collapse = "")
names(Genome) <- names(GenomeList)
SampleVCF <- readVcf(file = "D:/L1polymORF/Data/chr3_1161_aln.vcf", 
                     genome = Genome)
VCFDF <- as.data.frame(SampleVCF@fixed)

gVCF <- geno(SampleVCF)

PLDF <- as.data.frame(gVCF$PL)

sapply(geno(SampleVCF), class)
snpSummary(SampleVCF)

as.data.frame(info(header(SampleVCF)))
info(SampleVCF)
head(rowRanges(SampleVCF))

pass(SampleVCF)
any(qual(SampleVCF) != 0)
qual
