# This script prepares SNP data for NA 12878 to get snps overlapping
# with regions sequenced by PacBio capture data

# Load packages
library(seqinr)
library(ShortRead)
library(Rsamtools)
library(rtracklayer)
library(csaw)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Minimum proportion of SNP in reads to be called
MinPropCall <- 0.6

# Paths (SCG4)
PathStart             <- '/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R'
PathChrLen            <- '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/ChromLengthsHg19.Rdata'
PathReadRanges_HiFi   <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioHiFiSubreads_L1Ranges.RData"
PathReadRanges_Normal <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_L1Ranges.RData"
PathSNP_NA12878       <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf"
PathBam_HiFi          <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/PacBioL1captPS/PacBioHiFiSubreads_hg19masked_ngmlr.sorted.bam"
PathBam_Normal        <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/BZ_NA12878L1capt5-9kb_subreads_hg19masked.sorted.bam"
PathOutput            <- "/srv/gsfs0/projects/levinson/hzudohna/PacBioCapture/SNPinRangesNA12878.RData"

# Source start script
source(PathStart)

###############################################
#                                             #
#    Load ranges where reads were mapped      #
#                                             #
###############################################

# Load GRanges from HiFi polymerase data
load(PathReadRanges_HiFi)
IslGRanges_reduced_HiFi <- IslGRanges_reduced

# Load GRanges from HiFi polymerase data
load(PathReadRanges_Normal)
IslGRanges_reduced_Normal <- IslGRanges_reduced

# Get intersection of both sets of genomic range
GRUnion <- union(IslGRanges_reduced_HiFi, IslGRanges_reduced_Normal)

# Load data on chromosome length
load(PathChrLen)

###############################################
#                                             #
#    Get SNPs overlapping with above ranges   #
#                                             #
###############################################

# Get SNPs of NA12878 overlapping with genomic ranges
VCFLines <- readLines(PathSNP_NA12878, n = 200)

SNP_GR       <- GRanges()
NrLines2Read <- 10^5
LinesSkip    <- max(grep("##", VCFLines)) + 1
NrLinesRead  <- NrLines2Read
cat("*******   Subsetting SNPs   ********\n")
while(NrLinesRead >= NrLines2Read){
  SNPs <- read.delim(file = PathSNP_NA12878, skip = LinesSkip, nrows = NrLines2Read,
                    header = F)
  SNPs <- SNPs[,c(1:2, 2, 4, 5)]
  colnames(SNPs)  <- c("chromosome", "start", "end", "REF", "ALT")
  SNPs$chromosome <- paste("chr", SNPs$chromosome, sep = "")
  SNP_GR_local    <- makeGRangesFromDataFrame(SNPs, keep.extra.columns = T)
  SNP_GR_local    <- subsetByOverlaps(SNP_GR_local, GRUnion)
  SNP_GR          <- c(SNP_GR, SNP_GR_local)
  NrLinesRead     <- nrow(SNPs)
  LinesSkip       <- LinesSkip + NrLinesRead
  cat("Processed", LinesSkip, "lines\n")
}
rm(list = "SNPs")
gc()

###############################################
#                                             #
#       Get reads for above ranges            #
#                                             #
###############################################

cat("*******   Getting reads from two datasets   ********\n")

# Subset GRUnion to get the ones containing SNPs 
GRUnion_withSNP <- subsetByOverlaps(GRUnion, SNP_GR)

# Get reference sequence for each range
RefSeqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRUnion_withSNP)

# Get reads for different datasets
scanPars        <- ScanBamParam(which = GRUnion_withSNP, what = scanBamWhat())
ReadList_HiFi   <- scanBam(PathBam_HiFi, param = scanPars)
ReadList_Normal <- scanBam(PathBam_Normal, param = scanPars)

# Saving results
cat("*****  Saving results to", PathOutput, " ******\n")
save.image(PathOutput)