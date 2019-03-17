# The following script compares a vcf file of identified L1 insertions from
# a simulated genome with the insertions generated in the simulation

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
VcfDirs <- list.dirs(SimDir, full.names = F)
VcfDirs <- VcfDirs[VcfDirs %in% SampleColumns]
VcfDirs <-  setdiff(VcfDirs, c("HG00096", "HG00097", "HG00099", "HG00100","HG00101"))
VcfFiles <- paste(SimDir, VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# VcfDirs <-  c("HG00096", "HG00097", "HG00099", "HG00100","HG00101")
# VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
#       VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# Read in vcf file and create genomic ranges
L1Detect <- NULL
for (VcfFile in VcfFiles[file.exists(VcfFiles)]){
  cat("Processing", VcfFile, "\n")
  
  VcfFile <- ReadVCF(VcfFile)
  VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                      start.field = "POS",
                                      end.field = "POS")
  VcfGR <- resize(VcfGR, width = 150, fix = "center")
  
  # Get sample ID from vcf column and get index of L1 insertions that occur in
  # that sample
  SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
  cat("Processing", SampleID, "\n")
  idxL1    <- which(L1_1000G[,SampleID] > 0)
  SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
                                       start.field = "POS",
                                       end.field = "POS")
  
  # Create data.frame that keeps track of L1 width and their detection status
  L1DetectNew <- data.frame(blnDetect = overlapsAny(SampleGR, VcfGR),
                            L1width = L1_1000G$InsLength[idxL1])
  L1Detect   <- rbind(L1Detect, L1DetectNew)
  
}
  
# Add column for full-length L1
L1Detect$blnFull <- L1Detect$L1width >= 6000

# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width <- glm(blnDetect ~ L1width + blnFull, data = L1Detect,
                            family = binomial)
summary(LogReg_DetectL1width)


# Plot the percentage of detected LINE-1 per length class
# Create a vector of L1 start classes
L1Detect$InsLengthClass <- cut(L1Detect$L1width, breaks = 
                                  seq(0, 6500, 500))

# Get mean L1 detection per L1 width
L1WidthAggregated <- aggregate(L1Detect[,c("L1width", "blnDetect")], 
                               by = list(L1Detect$InsLengthClass), 
                               FUN = function(x) mean(x, na.rm = T))
L1WidthAggregated_n <- aggregate(L1Detect$L1width, 
                                 by = list(L1Detect$InsLengthClass), 
                                 FUN = function(x) sum(!is.na(x)))
L1WidthAgg <- merge(L1WidthAggregated, L1WidthAggregated_n)

# Plot detection probability vs insertion size
# pdf(file = "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/DetectVsInsertSize.pdf")
# plot(L1WidthAgg$L1width, L1WidthAgg$blnDetect, xlab = "L1 width",
#      ylab = "Proportion detected")
# SDdetect <- sqrt(L1WidthAgg$blnDetect*(1 - L1WidthAgg$blnDetect) /
#   L1WidthAgg$x)
# AddErrorBars(MidX = L1WidthAgg$L1width, 
#              MidY = L1WidthAgg$blnDetect, 
#              ErrorRange = SDdetect,
#              TipWidth = 20)
# dev.off()

#############################################
#                                           #
#   Analyze L1 detection from filtered      #
#          simulated files                  #
#                                           #
#############################################

# Get names of vcf files
VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/",
                     full.names = F)
VcfDirs <-  VcfDirs[VcfDirs %in% paste(SampleColumns, "Filtered", sep = "_")]
VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
                  VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# VcfDirs <-  c("HG00096", "HG00097", "HG00099", "HG00100","HG00101")
# VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
#       VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# Read in vcf file and create genomic ranges
VcfFiles[file.exists(VcfFiles)]
L1Detect_Filter <- NULL
for (VcfFile in VcfFiles[file.exists(VcfFiles)]){
  cat("Processing", VcfFile, "\n")
  
  if(file.info(VcfFile)$size > 0){
    cat("Appending L1 calls\n")
    VcfFile <- ReadVCF(VcfFile)
    VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                        start.field = "POS",
                                        end.field = "POS")
    VcfGR <- resize(VcfGR, width = 150, fix = "center")
    
    # Get sample ID from vcf column and get index of L1 insertions that occur in
    # that sample
    SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
    cat("Processing", SampleID, "\n")
    idxL1    <- which(L1_1000G[,SampleID] > 0)
    SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
                                         start.field = "POS",
                                         end.field = "POS")
    
    # Create data.frame that keeps track of L1 width and their detection status
    L1DetectNew <- data.frame(blnDetect = overlapsAny(SampleGR, VcfGR),
                              L1width = L1_1000G$InsLength[idxL1])
    L1Detect_Filter   <- rbind(L1Detect_Filter, L1DetectNew)
    
  } else {
    cat("Empty file!\n")
  }
  
}


# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width_Filter <- glm(blnDetect ~ L1width, data = L1Detect_Filter,
                            family = binomial)
summary(LogReg_DetectL1width_Filter)


######################################################
#                                                    #
#   Analyze L1 detection from group analysis of      #
#          simulated files                           #
#                                                    #
######################################################

VcfFile <- ReadVCF(paste(SimDir, "LINE1.final_comp.vcf", sep = ""))
VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                    start.field = "POS",
                                    end.field = "POS")
VcfGR <- resize(VcfGR, width = 150, fix = "center")

# Get sample ID from vcf column and get index of L1 insertions that occur in
# that sample
SampleIDs <- sapply(10:ncol(VcfFile), function(x) {
  strsplit(colnames(VcfFile)[x], "_")[[1]][2]
})
L1Detect_Group <- NULL
for (SampleID in  SampleIDs){
  cat("Processing", SampleID, "\n")
  idxL1    <- which(L1_1000G[,SampleID] > 0)
  SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
                                       start.field = "POS",
                                       end.field = "POS")
  
  # Create data.frame that keeps track of L1 width and their detection status
  L1DetectNew <- data.frame(blnDetect = overlapsAny(SampleGR, VcfGR),
                            L1width = L1_1000G$InsLength[idxL1])
  L1Detect_Group <- rbind(L1Detect_Group, L1DetectNew)
  
}

# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width_Group <- glm(blnDetect ~ L1width, data = L1Detect_Group,
                            family = binomial)
summary(LogReg_DetectL1width_Group)


######################################################
#                                                    #
#             Save objects                           #
#                                                    #
######################################################

save.image(file = "/labs/dflev/hzudohna/1000Genomes/L1Simulated_MELT.RData")
