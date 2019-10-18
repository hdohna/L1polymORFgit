# The script below estimates selection coefficients of L1 from the 
# 1000 genome data

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

library(pracma)
# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load 1000 genome data
load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")

#####################################################################
#                                                                   #
#     Plot theoretical probability densities of frequencies         #
#                                                                   #
#####################################################################

# Create vectors of frequency values, selection coefficients, and colors
FreqVals <- seq(0, 1, 0.001)
SCoeffs  <- seq(-10^-3, 10^-3, 5*10^-4)
Cols     <- rainbow(length(SCoeffs))

jpeg(file = "D:/OneDrive - American University of Beirut/Teaching/Population genetics Fall 2019/Slides/Scoeff.jpeg",
     width = 700, height = 500, quality = 100)

# Create values of probability density (Distn1) and plot then
Distn1   <- sapply(FreqVals, function(x) AlleleFreqDistn(y = x, s = 0, N = 10^5))
plot(FreqVals, Distn1, type = "l", ylim = c(0, 20), xlim = c(0, 1),
     xlab = "Allele frequency",
     ylab = "Probability density")
for (i in 1:length(SCoeffs)){
  Distn <- sapply(FreqVals, function(x) AlleleFreqDistn(y = x, s = SCoeffs[i], N = 10^5))
  lines(FreqVals, Distn, col = Cols[i])
}
legend("topright", legend = paste("s =", SCoeffs), col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.9)
dev.off()
CreateDisplayPdf('D:/L1polymORF/Figures/FreqDistPerSelectCoeff.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

#####################################################################
#                                                                   #
#       Plot distribution of observed frequencies in a sample       #
#                                                                   #
#####################################################################

CountVals <- 1:5000
DistnCount1   <- sapply(CountVals, function(x) AlleleFreqSample(x, 
                  s = SCoeffs[1], N = 10^5, SampleSize = 5001, 
                  DetectProb = 0.8,
                  LogRegCoeff = c(-4, 9), blnIns = T))

plot(CountVals, exp(DistnCount1), type = "l", 
     xlab = "Allele frequency",
     ylab = "Probability density", xlim = c(0, 100))

for (i in 1:length(SCoeffs)){
  Distn <- sapply(CountVals, function(x) AlleleFreqSample(x, 
                 s = SCoeffs[i], N = 10^5, SampleSize = 5001, 
                 DetectProb = 0.8, LogRegCoeff = c(-4, 9), blnIns = T))
  lines(CountVals, exp(Distn), col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)

#####################################################################
#                                                                   #
#       Plot distribution of observed frequencies in a sample       #
#       simplified version 1: no accounting of detection
#                                                                   #
#####################################################################

DistnCount1   <- sapply(CountVals, function(x) AlleleFreqSample_simplified1(x, 
                    s = SCoeffs[1], N = 10^5, SampleSize = 5001, 
                    DetectProb = 0.8,
                    LogRegCoeff = c(-4, 9), blnIns = T))

plot(CountVals, exp(DistnCount1), type = "l", 
     xlab = "Allele frequency",
     ylab = "Probability density", xlim = c(0, 100))

for (i in 1:length(SCoeffs)){
  Distn <- sapply(CountVals, function(x) AlleleFreqSample_simplified1(x, 
                       s = SCoeffs[i], N = 10^5, SampleSize = 5001, 
                       DetectProb = 0.8, LogRegCoeff = c(-4, 9), blnIns = T))
  lines(CountVals, exp(Distn), col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)


#####################################################################
#                                                                   #
#       Plot distribution of observed frequencies in a sample       #
#       simplified version 2: no accounting of insertion
#                                                                   #
#####################################################################

DistnCount1   <- sapply(CountVals, function(x) AlleleFreqSample_simplified2(x, 
                        s = SCoeffs[1], N = 10^5, SampleSize = 5001, 
                       DetectProb = 0.8,
                       LogRegCoeff = c(-4, 9), blnIns = T))

plot(CountVals, exp(DistnCount1), type = "l", 
     xlab = "Allele frequency",
     ylab = "Probability density", xlim = c(0, 100))

for (i in 1:length(SCoeffs)){
  Distn <- sapply(CountVals, function(x) AlleleFreqSample_simplified2(x, 
                                                                      s = SCoeffs[i], N = 10^5, SampleSize = 5001, 
                                                                      DetectProb = 0.8, LogRegCoeff = c(-4, 9), blnIns = T))
  lines(CountVals, exp(Distn), col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)


#####################################################################
#                                                                   #
#       Plot distribution of observed frequencies in a sample       #
#       simplified version 3: no accounting of insertion
#       or detection
#                                                                   #
#####################################################################

SCoeffs  <- seq(-10^-4, 10^-4, 5*10^-5)
SCoeffs  <- seq(-10^-3, 10^-3, 5*10^-4)
SCoeffs  <- seq(-10^-2, 10^-2, 5*10^-3)
SCoeffs  <- seq(-5*10^-1, 5*10^-1, 10^-1)
Cols     <- rainbow(length(SCoeffs))
DistnCount1   <- sapply(CountVals, function(x) AlleleFreqSample_simplified3(x, 
                                   s = SCoeffs[1], N = 10^6, SampleSize = 5001, 
                                   DetectProb = 0.5,
                                   LogRegCoeff = c(-4, 9), blnIns = T))

plot(CountVals, DistnCount1, type = "l", 
     xlab = "Allele frequency",
     ylab = "Probability density", xlim = c(0, 100))

for (i in 1:length(SCoeffs)){
  Distn <- sapply(CountVals, function(x) AlleleFreqSample_simplified3(x, 
               s = SCoeffs[i], N = 10^6, SampleSize = 5001, 
              DetectProb = 0.5, LogRegCoeff = c(-4, 9), blnIns = T))
  lines(CountVals, Distn, col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)

#####################################################################
#                                                                   #
#       Plot distribution of observed frequencies in a sample       #
#       simplified version 4: no accounting of insertion
#       or detection and different numeric integration
#                                                                   #
#####################################################################

SCoeffs  <- seq(-10^-4, 10^-4, 5*10^-5)
SCoeffs  <- seq(-10^-3, 10^-3, 5*10^-4)
SCoeffs  <- seq(-10^-2, 10^-2, 5*10^-3)
# SCoeffs  <- seq(-5*10^-1, 5*10^-1, 10^-1)
Cols     <- rainbow(length(SCoeffs))
DistnCount1   <- sapply(CountVals, function(x) AlleleFreqSample_simplified4(x, 
                                                                            s = SCoeffs[1], N = 10^6, SampleSize = 5001, 
                                                                            DetectProb = 0.5,
                                                                            LogRegCoeff = c(-4, 9), blnIns = T))

plot(CountVals, DistnCount1, type = "l", 
     xlab = "Allele frequency",
     ylab = "Probability density", xlim = c(0, 100))

for (i in 1:length(SCoeffs)){
  Distn <- sapply(CountVals, function(x) AlleleFreqSample_simplified4(x, 
                                                                      s = SCoeffs[i], N = 10^6, SampleSize = 5001, 
                                                                      DetectProb = 0.5, LogRegCoeff = c(-4, 9), blnIns = T))
  lines(CountVals, Distn, col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)

#####################################################################
#                                                                   #
#       Histograms of actual counts
#                                                                   #
#####################################################################

# 1000 genome data
L1_1000G$FreqCount <- rowSums(L1_1000G[,SampleColumns])
hist(L1_1000G$FreqCount, breaks = 0:5000, xlim = c(0, 50))
hist(L1_1000G$FreqCount[L1_1000G$InsLength <= 500], breaks = 0:5000, xlim = c(0, 50))

# Read in vcf file with MELT insertion calls
MEInsCall <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf", 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEInsCall <- MEInsCall[MEInsCall$Type == "<INS:ME:LINE1>",]

# Extract allele frequency from info column
GetAF <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  AFch   <- strsplit(xSplit[length(xSplit)], "=")[[1]][2]
  as.numeric(AFch)
}
GetLength <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  LengthCh   <- strsplit(xSplit[grep("SVLEN=", xSplit)], "=")[[1]][2]
  as.numeric(LengthCh)
}

# Add columns necessary for analysis 
MEInsCall$AF <- sapply(MEInsCall$Info, GetAF)
MEInsCall <- MEInsCall[!is.na(MEInsCall$AF), ]
MEInsCall$L1width <- sapply(MEInsCall$Info, GetLength)
# MEInsCall$SampleSize <- 2 * MEInsSamplesize
MEInsCall$Freq <- ceiling(MEInsCall$SampleSize * MEInsCall$AF) # TODO: Figure out why not integers!
MEInsCall$blnFull <- MEInsCall$L1width >= MinLengthFullL1

hist(MEInsCall$Freq, breaks = 0:5000, xlim = c(0, 50))
hist(MEInsCall$Freq[MEInsCall$L1width <= 500], breaks = 0:5000, xlim = c(0, 50))
