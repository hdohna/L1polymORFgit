##############################################
#
# General description:
#
#   The following script reads bam files from my read mapping and bed files from
#   Bo's mapping and peak find and compares them

# Input:
#
#    SRIP_eul1db: data on methods

# Output:
#   
#    : ...

##############################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(ShortRead)
library(csaw)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg38)

#######################################
#                                     #
#    Turn BAM files into GRanges      #
#                                     #
#######################################
list.files("D:/L1polymORF/Data/")
cat("*******   Turning BAM files into GRanges ...   *******\n")

InputFile = "D:/L1polymORF/Data/NA12878-L1HS_S1_L001_R1_001_16b02f092319.bam"
Chrom = "chr1"
MinPerIsland = 1
ChromLength <- length(BSgenome.Hsapiens.UCSC.hg38[[Chrom]])
R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength))

Reads <- extractReads(bam.file = InputFile, region = R1, 
                      param=readParam(pe = "both", max.frag = 10000))
ReadCov <- coverage(Reads)
Islands <- slice(ReadCov, lower = MinPerIsland)

cat("*******   Turning BED files into GRanges ...   *******\n")
L1Ref        <- import.bed(con = "D:/L1polymORF/Data/hg38.fa_L1HS_6kb.bed") 
NA12878Peaks <- import.bed(con = "D:/L1polymORF/Data/NA12878_L1_capt_IDTx_peaks.bed") 

#
NA12878Peaks_chr1 <- NA12878Peaks[seqnames(NA12878Peaks) == "chr1"]
sum(NA12878Peaks_chr1 %over% Islands) / length(NA12878Peaks_chr1)



