library("VariantAnnotation")
library("seqinr")


GenomeList <- read.fasta("D:/L1polymORF/Data/Homo_sapiens_L1_consensus.fa")
Genome <- paste(GenomeList[[1]], collapse = "")
names(Genome) <- names(GenomeList)
SampleVCF <- readVcf(file = "D:/L1polymORF/Data/chr3_1161_aln.vcf", 
                     genome = Genome)
VCFDF <- as.data.frame(SampleVCF@fixed)
geno(SampleVCF)
snpSummary(SampleVCF)
