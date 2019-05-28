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
#SimDir <- "D:/L1polymORF/Data/SimulatedL1/"

# Get names of vcf files
VcfDirs <- list.dirs(SimDir, full.names = F)
VcfDirs <- VcfDirs[VcfDirs %in% SampleColumns]
#VcfDirs <-  setdiff(VcfDirs, c("HG00096", "HG00097", "HG00099", "HG00100","HG00101"))
VcfFiles <- paste(SimDir, VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# VcfDirs <-  c("HG00096", "HG00097", "HG00099", "HG00100","HG00101")
# VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
#       VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# Read in vcf file and create genomic ranges
L1Detect <- NULL
VcfFile <- VcfFiles[1]
for (VcfFile in VcfFiles[file.exists(VcfFiles)]){
  cat("Processing", VcfFile, "\n")
  
  # Read in vcf file and make genomic ranges
  VcfFile <- ReadVCF(VcfFile)
  VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                      start.field = "POS",
                                      end.field = "POS")
  VcfGR <- resize(VcfGR, width = 50, fix = "center")
  
  # Get estimated L1 length and genotype from vcf file
  L1widthVcf <- sapply(VcfFile$INFO, GetFromVcfINFO_SVLength)
  L1GenoVcf  <- sapply(VcfFile[,ncol(VcfFile)], GetFromVcfGeno_GenoNum)
  L1StartEndVcf  <- t(sapply(VcfFile[,ncol(VcfFile)], GetFromVcfINFO_MELT_L1StartEnd))
  
  # Get sample ID from vcf column and get index of L1 insertions that occur in
  # that sample
  SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
  cat("Processing", SampleID, "\n")
  idxL1    <- which(L1_1000G[,SampleID] > 0)
  SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
                                       start.field = "POS",
                                       end.field = "POS")
  
  # Determine overlaps between simulated and detected L1
  OLSimDetect <- findOverlaps(VcfGR, SampleGR)
  
  # Create data.frame that keeps track of L1 present in 1000 Genomes and 
  # their detection status
  L1DetectNew <- data.frame(Chrom  = L1_1000G$CHROM[idxL1],
                            PosTrue  = L1_1000G$POS[idxL1],
                            PosEst      = NA,
                            blnDetect   = overlapsAny(SampleGR, VcfGR),
                            L1GenoTrue  = L1_1000G[idxL1, SampleID],
                            L1widthTrue = L1_1000G$InsLength[idxL1],
                            L1StartTrue = L1_1000G$L1Start[idxL1],
                            L1EndTrue   = L1_1000G$L1End[idxL1],
                            L1widthEst  = NA,
                            L1StartEst  = NA,
                            L1EndEst  = NA,
                            L1GenoEst   = NA,
                            EstFilter   = NA,
                            SampleID = SampleID)
  L1DetectNew$PosEst[OLSimDetect@to]     <- VcfFile$POS[OLSimDetect@from]
  L1DetectNew$L1widthEst[OLSimDetect@to] <- L1widthVcf[OLSimDetect@from]
  L1DetectNew$L1StartEst[OLSimDetect@to] <- L1StartEndVcf[OLSimDetect@from, 1]
  L1DetectNew$L1EndEst[OLSimDetect@to]   <- L1StartEndVcf[OLSimDetect@from, 2]
  L1DetectNew$L1GenoEst[OLSimDetect@to]  <- L1GenoVcf[OLSimDetect@from]
  L1DetectNew$EstFilter[OLSimDetect@to]  <- VcfFile$FILTER[OLSimDetect@from]
  L1Detect   <- rbind(L1Detect, L1DetectNew)
  
  # Create data.frame that keeps track of L1 not present in 1000 Genomes and 
  # their detection status
  blnPresent  <- overlapsAny(VcfGR, SampleGR)
  if(any(!blnPresent)){
    L1DetectNew <- data.frame(Chrom  = VcfFile$X.CHROM[!blnPresent],
                              PosTrue  = NA,
                              PosEst      = VcfFile$POS[!blnPresent],
                              blnDetect   = NA,
                              L1GenoTrue  = 0,
                              L1widthTrue = 0,
                              L1StartTrue = NA,
                              L1EndTrue   = NA,
                              L1widthEst  = L1widthVcf[!blnPresent],
                              L1StartEst  =  L1StartEndVcf[!blnPresent, 1],
                              L1EndEst    =  L1StartEndVcf[!blnPresent, 2],
                              L1GenoEst   = L1GenoVcf[!blnPresent],
                              EstFilter   = VcfFile$FILTER[!blnPresent],
                              SampleID = SampleID)
    L1Detect   <- rbind(L1Detect, L1DetectNew)
  }
}
  
# Add column for full-length L1
L1Detect$blnFull <- L1Detect$L1widthTrue >= 6000

# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width <- glm(blnDetect ~ L1widthTrue + blnFull, data = L1Detect,
                            family = binomial)
summary(LogReg_DetectL1width)

# Return some summaries
Sensitivity <- mean(L1Detect$blnDetect, na.rm = T)
blnL1Estimated <- !is.na(L1Detect$L1GenoEst)
Specificity   <- sum(blnL1Estimated & L1Detect$blnDetect, na.rm = T) / 
  sum(blnL1Estimated)
cat("Sensitivity:", Sensitivity, "\n")
cat("Specificity:", Specificity, "\n")

# Plot estimated vs true L1 length
plot(L1Detect$L1widthTrue, L1Detect$L1widthEst)
cor.test(L1Detect$L1widthTrue, L1Detect$L1widthEst, method = "spearman")
FullConTab <- table(L1Detect$blnFull, L1Detect$L1widthEst >= 6000)
t(t(FullConTab) / colSums(FullConTab))
chisq.test(L1Detect$blnFull, L1Detect$L1widthEst >= 6000)
hist(L1Detect$L1widthEst)

# Probability of false positive, given the estimated genotype (0, 1, or 2)
GenoTrueEst <- table(L1Detect$L1GenoTrue, L1Detect$L1GenoEst)
GenoTrueEst[1,] / colSums(GenoTrueEst)
barplot(GenoTrueEst[1,] / colSums(GenoTrueEst), xlab = "Estimated genotype",
        ylab = )

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

###################################################
#                                                 #
#   Analyze L1 detection from simulated L1        #
#     with variable difference from consensus     #
#                                                 #
###################################################

# # Load file with simulation specifications
# load("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/L1SimSettings.RData")
# 
# # Specify simulation directory
# SimDir <- "/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/"
# 
# # Get names of vcf files
# VcfDirs  <- list.dirs(SimDir, full.names = F)
# VcfDirs  <- VcfDirs[VcfDirs %in% paste(SampleColumns, "Var", sep = "_")]
# VcfFiles <- paste(SimDir, VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# L1DetectVar <- NULL
# for (VcfFile in VcfFiles[file.exists(VcfFiles)]){
#   cat("Processing", VcfFile, "\n")
#   
#   VcfFile <- ReadVCF(VcfFile)
#   VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
#                                       start.field = "POS",
#                                       end.field = "POS")
#   VcfGR <- resize(VcfGR, width = 150, fix = "center")
#   
#   # Get sample ID from vcf column and get index of L1 insertions that occur in
#   # that sample
#   SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
#   cat("Processing", SampleID, "\n")
#   idxL1    <- which(L1_1000G[,SampleID] > 0)
#   SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
#                                        start.field = "POS",
#                                        end.field = "POS")
#   
#   # Create data.frame that keeps track of L1 width and their detection status
#   L1DetectNew <- data.frame(blnDetect = overlapsAny(SampleGR, VcfGR),
#                             L1PropDiff = L1_1000G$PropDiff[idxL1])
#   L1DetectVar   <- rbind(L1DetectVar, L1DetectNew)
#   
# }
# 
# 
# # Logistic regression for detection probability as finction of insertion length
# LogReg_DetectL1L1PropDiff <- glm(blnDetect ~ L1PropDiff, data = L1Detect,
#                             family = binomial)
# summary(LogReg_DetectL1L1PropDiff)
# 
# 
# # Plot the percentage of detected LINE-1 per length class
# # Create a vector of L1 start classes
# L1Detect$InsLengthClass <- cut(L1Detect$L1width, breaks = 
#                                  seq(0, 6500, 500))
# 
# # Get mean L1 detection per L1 width
# L1WidthAggregated <- aggregate(L1Detect[,c("L1width", "blnDetect")], 
#                                by = list(L1Detect$InsLengthClass), 
#                                FUN = function(x) mean(x, na.rm = T))
# L1WidthAggregated_n <- aggregate(L1Detect$L1width, 
#                                  by = list(L1Detect$InsLengthClass), 
#                                  FUN = function(x) sum(!is.na(x)))
# L1WidthAgg <- merge(L1WidthAggregated, L1WidthAggregated_n)
# 
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

# # Get names of vcf files
# VcfDirs <- list.dirs("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/",
#                      full.names = F)
# VcfDirs <-  VcfDirs[VcfDirs %in% paste(SampleColumns, "Filtered", sep = "_")]
# VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
#                   VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# # VcfDirs <-  c("HG00096", "HG00097", "HG00099", "HG00100","HG00101")
# # VcfFiles <- paste("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/", 
# #       VcfDirs, "/LINE1.final_comp.vcf", sep = "")
# # Read in vcf file and create genomic ranges
# VcfFiles[file.exists(VcfFiles)]
# L1Detect_Filter <- NULL
# for (VcfFile in VcfFiles[file.exists(VcfFiles)]){
#   cat("Processing", VcfFile, "\n")
#   
#   if(file.info(VcfFile)$size > 0){
#     cat("Appending L1 calls\n")
#     VcfFile <- ReadVCF(VcfFile)
#     VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
#                                         start.field = "POS",
#                                         end.field = "POS")
#     VcfGR <- resize(VcfGR, width = 150, fix = "center")
#     
#     # Get sample ID from vcf column and get index of L1 insertions that occur in
#     # that sample
#     SampleID <- strsplit(colnames(VcfFile)[ncol(VcfFile)], "_")[[1]][2]
#     cat("Processing", SampleID, "\n")
#     idxL1    <- which(L1_1000G[,SampleID] > 0)
#     SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], seqnames.field = "chromosome",
#                                          start.field = "POS",
#                                          end.field = "POS")
#     
#     # Create data.frame that keeps track of L1 width and their detection status
#     L1DetectNew <- data.frame(blnDetect = overlapsAny(SampleGR, VcfGR),
#                               L1width = L1_1000G$InsLength[idxL1])
#     L1Detect_Filter   <- rbind(L1Detect_Filter, L1DetectNew)
#     
#   } else {
#     cat("Empty file!\n")
#   }
#   
# }
# 
# 
# # Logistic regression for detection probability as finction of insertion length
# LogReg_DetectL1width_Filter <- glm(blnDetect ~ L1width, data = L1Detect_Filter,
#                             family = binomial)
# summary(LogReg_DetectL1width_Filter)
# 

######################################################
#                                                    #
#   Analyze L1 detection from group analysis of      #
#          simulated files                           #
#                                                    #
######################################################

VcfFile <- ReadVCF("/labs/dflev/hzudohna/1000Genomes/L1_simulation_MELT/LINE1.final_comp.vcf")
VcfGR   <- makeGRangesFromDataFrame(VcfFile, seqnames.field = "X.CHROM",
                                    start.field = "POS",
                                    end.field = "POS")
VcfGR <- resize(VcfGR, width = 100, fix = "center")

# Get estimated L1 length and genotype from vcf file
L1widthVcf <- sapply(VcfFile$INFO, GetFromVcfINFO_SVLength)
L1GenoVcf  <- sapply(VcfFile[,ncol(VcfFile)], GetFromVcfGeno_GenoNum)
L1StartEndVcf  <- t(sapply(VcfFile[,ncol(VcfFile)], GetFromVcfINFO_MELT_L1StartEnd))

# Get sample ID from vcf column and get index of L1 insertions that occur in
# that sample
SampleIDs <- sapply(10:ncol(VcfFile), function(x) {
  strsplit(colnames(VcfFile)[x], "_")[[1]][2]
})
L1Detect_Group <- NULL
for (SampleID in  SampleIDs){
  cat("Processing", SampleID, "\n")
  idxL1    <- which(L1_1000G[,SampleID] > 0)
  SampleGR <- makeGRangesFromDataFrame(L1_1000G[idxL1,], 
                                       seqnames.field = "chromosome",
                                       start.field = "POS",
                                       end.field = "POS")
  
  # Determine overlaps between simulated and detected L1
  OLSimDetect <- findOverlaps(VcfGR, SampleGR)
  
  # Create data.frame that keeps track of L1 width and their detection status
  L1DetectNew <- data.frame(Chrom  = L1_1000G$CHROM[idxL1],
                            PosTrue  = L1_1000G$POS[idxL1],
                            PosEst   = NA,
                            blnDetect   = overlapsAny(SampleGR, VcfGR),
                            L1GenoTrue  = L1_1000G[idxL1, SampleID],
                            L1widthTrue = L1_1000G$InsLength[idxL1],
                            L1StartTrue = L1_1000G$L1Start[idxL1],
                            L1EndTrue   = L1_1000G$L1End[idxL1],
                            L1widthEst  = NA,
                            L1StartEst  = NA,
                            L1EndEst  = NA,
                            L1GenoEst   = NA,
                            EstFilter   = NA,
                            SampleID = SampleID)
  L1DetectNew$PosEst[OLSimDetect@to]     <- VcfFile$POS[OLSimDetect@from]
  L1DetectNew$L1widthEst[OLSimDetect@to] <- L1widthVcf[OLSimDetect@from]
  L1DetectNew$L1StartEst[OLSimDetect@to] <- L1StartEndVcf[OLSimDetect@from, 1]
  L1DetectNew$L1EndEst[OLSimDetect@to]   <- L1StartEndVcf[OLSimDetect@from, 2]
  L1DetectNew$L1GenoEst[OLSimDetect@to]  <- L1GenoVcf[OLSimDetect@from]
  L1DetectNew$EstFilter[OLSimDetect@to]  <- VcfFile$FILTER[OLSimDetect@from]
  L1Detect_Group <- rbind(L1Detect_Group, L1DetectNew)
  
  # Create data.frame that keeps track of L1 not present in 1000 Genomes and 
  # their detection status
  blnPresent  <- overlapsAny(VcfGR, SampleGR)
  if(any(!blnPresent)){
    L1DetectNew <- data.frame(Chrom  = VcfFile$X.CHROM[!blnPresent],
                              PosTrue  = NA,
                              PosEst      = VcfFile$POS[!blnPresent],
                              blnDetect   = NA,
                              L1GenoTrue  = 0,
                              L1widthTrue = 0,
                              L1StartTrue = NA,
                              L1EndTrue   = NA,
                              L1widthEst  = L1widthVcf[!blnPresent],
                              L1StartEst  =  L1StartEndVcf[!blnPresent, 1],
                              L1EndEst    =  L1StartEndVcf[!blnPresent, 2],
                              L1GenoEst   = L1GenoVcf[!blnPresent],
                              EstFilter   = VcfFile$FILTER[!blnPresent],
                              SampleID = SampleID)
  }

  L1Detect_Group <- rbind(L1Detect_Group, L1DetectNew)
  
}

# Logistic regression for detection probability as finction of insertion length
LogReg_DetectL1width_Group <- glm(blnDetect ~ L1widthTrue, data = L1Detect_Group,
                            family = binomial)
summary(LogReg_DetectL1width_Group)


######################################################
#                                                    #
#             Save objects                           #
#                                                    #
######################################################

save.image(file = "/labs/dflev/hzudohna/1000Genomes/L1Simulated_MELT.RData")
