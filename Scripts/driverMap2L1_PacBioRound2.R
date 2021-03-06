##############################################
#
# General description:
#
#   The following script reads in ranges of non-reference L1 peaks, and a 
#   PacBio bam file, writes out a fastq file aligns all fastq files in a folder
#   to the same reference using bwa

# Input:
#
#     OutFastQFolder: folder that contains the fastq files to be mapped
#     Homo_sapiens_L1_consensus.fa: L1 consensus sequence

# Output:
#   
#    Little sam file for each fastq file

##############################################

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Load necessary libraries
library(GenomicRanges)
library(Rsamtools)
library(csaw)

# Specify file paths 
BamFilePacBio     <- "/home/hzudohna/sorted_final_merged.bam"
OutFastQFolder    <- "/home/hzudohna/L1polymORF/Data/PacBioFastqPerSuspectPeak_FromPacBio/"
L1Consensus       <- "/home/hzudohna/L1polymORF/Data/Homo_sapiens_L1_consensus.fa"
CoverSummaryPlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverNonRef_Pacbio_Round2.pdf'
CoverComparePlot  <- '/home/hzudohna/L1polymORF/Figures/L1HSCoverEx_Pacbio_Round2.pdf'
OutResults        <- '/home/hzudohna/L1polymORF/Data/L1NonReference_Pacbio_Round2.Rdata'

# Boolean indicators for different actions
blnRun_getNonRefL1  <- TRUE
blnWriteFastq       <- TRUE
blnMap2L1           <- TRUE
blnAddReadGroups    <- TRUE
blnCallHaplotypes   <- TRUE
blnAnalyze         <- TRUE

# Suffices for alignment files created by BWA
SamSuffix <- "_aln.sam"
BamSuffix <- paste(substr(SamSuffix, 1, nchar(SamSuffix) - 4), ".bam", sep = "")

# Run or load results from script 'getNonRefL1Ranges_bix2.R
#load("D:/L1polymORF/Data/NonRefL1Ranges.Rdata")
if (blnRun_getNonRefL1) {
  cat("*******   Running script getNonRefL1Ranges_bix2.R   *******\n")
  source('/home/hzudohna/L1polymORF/Scripts/AnalyzePacBioL1Ranges_bix2.R')
} else {
  load("/home/hzudohna/L1polymORF/Data/AnalyzedPacBioL1Ranges.RData")
}

#######################################################
#                                                     #
#    Write fastq of suspected L1 not in reference     #
#                                                     #
#######################################################

if(blnWriteFastq){
  cat("*******   Writing little fastq files ...   *******\n")
    
  cat("Reading Pacbio reads per range \n")
  param <- ScanBamParam(which = SuspectL1Ranges, 
                        what = scanBamWhat())
  ScannedReads <- scanBam(file = BamFilePacBio, 
                          param = param)
  
  # creat a vector of file prefixes
  FilePrefix <- paste(seqnames(SuspectL1Ranges), idxSuspectL1Ranges, 
                      sep = "_")
  if (any(duplicated(FilePrefix))){
    stop("Duplicated fastq file names")
  }

  FilePaths <- t(sapply (1:length(ScannedReads), function(i){
    cat("Writing reads for peak", i, "of", length(ScannedReads), "\n")
    Reads <- ScannedReads[[i]]
    if (length(Reads$seq) > 0){
      WriteFastqAndSample(Reads, OutFastQFolder, FilePrefix[i])
      cat("Written", length(Reads$seq),"reads to fastq file ", FilePrefix[i], "\n")
    } else {
      cat("No read in range - no fastq file written out\n")
    }
  }))
}

#######################################
#                                     #
#     Map fastq file per range        #
#            to L1HS                  #
#                                     #
#######################################

if(blnMap2L1){

    FilePaths <- MapMultiFastq(FastQFolder = OutFastQFolder,
     AlignCommand = '/home/txw/bwa/bwa-0.7.12/bwa mem -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0',
     Reference = L1Consensus)
}

#######################################
#                                     #
#     Add read groups                 #
#                                     #
#######################################

if(blnAddReadGroups){
  
  FilePathsRG <- AddMultiReadGroups(FastQFolder = OutFastQFolder,
                                    PicardCmd   = "java -jar /home/txw/picard/picard-tools-1.131/picard.jar AddOrReplaceReadGroups",
                                    OptionLines = c("RGLB=lib1", "RGPL=illumina", "RGPU=unit1", "RGSM=20",
                                                    "SORT_ORDER=null", "CREATE_INDEX=TRUE", "VALIDATION_STRINGENCY=LENIENT"))
}

#######################################
#                                     #
#     Call haplotypes                 #
#                                     #
#######################################

if(blnCallHaplotypes){
  
  FilePathsVCF <- CallMultiVariants(BamFolder = OutFastQFolder,  
                                    HapTypeCallCmd = "java -jar /home/txw/GATK/GenomeAnalysisTK-2.1-11-g13c0244/GenomeAnalysisTK.jar -T HaplotypeCaller",
                                    RefSeqPath = L1Consensus,
                                    OptionLines = "--emitRefConfidence GVCF",
                                    BamSuffix = "withRG.bam",
                                    VCFSuffix = ".vcf") 
}
#######################################
#                                     #
#     Import reads mapped to L1       #
#                                     #
#######################################

if (blnAnalyze){
  cat("*******  Reading and analyzing mapped reads ...   *******\n")
  
  # Get all names of sam files created by BWA
  #SamFileNames <- FilePaths$SamPath
  
  # Get all names of sam files created by BWA
  SamFileNames <- list.files(OutFastQFolder, pattern = SamSuffix,
                             full.names = T)
  
  # Turn sam files into bam files
  for (fn in SamFileNames) {
    cat("Turning", fn, "into a bam file\n")
    asBam(fn, destination = substr(fn, 1, nchar(fn) - 4), overwrite = T)
  }
  
  # get names of newly created bam files
  FileNames <- list.files(OutFastQFolder, pattern = BamSuffix,
                          full.names = T)
  FileNames <- FileNames[-grep(".bam.", FileNames)]
  
  # Loop through file names and read in bam files of reads mapped to L1
  ScannedL1Ranges <- lapply(FileNames, function(x) scanBam(x))
  
  # Count the number of reads mapped
  NrMapped2L1 <- sapply(ScannedL1Ranges, function(x){
    sum(!is.na(x[[1]]$pos))
  })
  
  # Get aligned reads per peak
  R1 <- GRanges(seqnames = "L1HS_L1_Homo_sapiens", 
                ranges = IRanges(start = 1, end = 6000))
  ReadsPerL1 <- lapply(FileNames[NrMapped2L1 > 0], function(x) {
    Reads <- extractReads(x, R1)
  })
  
  # Calculate a coverage matrix
  CoverMat <- t(sapply(ReadsPerL1, function(x){
    Cov <- coverage(x)
    as.vector(Cov[[1]])
  }))
  
  # Get means and 95% quantiles 
  QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
  idxFw <- 1:ncol(CoverMat)
  idxRv <- ncol(CoverMat):1
  pdf(file = CoverSummaryPlot)
  plot(QuantileMat[2,], type = "n", ylim = c(0, max(QuantileMat)), 
       ylab = 'Coverage', xlab = "Genomic position")
  polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
          col = "grey", border = NA)
  lines(QuantileMat[2,], lwd = 1.2)
  dev.off()
  
  # # Determine range index from file name 
  pdf(file = CoverComparePlot)
  plot(CoverMat[1,], type = "s", xlab = "Position on L1",
       ylab = "Coverage", ylim = c(0, 100))
  Cols <- rainbow(nrow(CoverMat))
  for (i in 1:nrow(CoverMat)){
    lines(CoverMat[i,], type = "s", col = Cols[i])
    
  }
  dev.off()
  
  # Save results
  cat("*******  Saving results ...   *******\n")
  save(list = c("IslGRanges_reduced", "FileNames", "ScannedL1Ranges", "ReadsPerL1", 
                "CoverMat", "QuantileMat"), file = OutResults)
  
}

