# The following script analyzes the correlation between L1 frequency and 
# insertion length in data collected by Ewing et al. 2011

# Load packages
library(ggplot2)
source("D:/L1polymORFgit/Scripts/_Start_L1polymORF.R")

# Load data by Ewing et al.2011
L1Ewing2011_S4 <- read.csv("D:/L1polymORF/Data/Ewing2011Table_S4.csv")

# Calculate frequency
N <- 2 * (ncol(L1Ewing2011_S4) - 6)
L1Ewing2011_S4$Frequency <- rowSums(L1Ewing2011_S4[,7:ncol(L1Ewing2011_S4)])/ N

# Group insertion length in units of 500 bp
InsLengthGroups <- seq(0, 6500, 500)
InsLengthMid <- InsLengthGroups[-length(InsLengthGroups)] + 250
L1Ewing2011_S4$InsLengthGroup <- cut(L1Ewing2011_S4$length,
                                       breaks = InsLengthGroups)
FreqMean <- aggregate(Frequency ~ InsLengthGroup, data = L1Ewing2011_S4, FUN = mean)
FreqMin <- aggregate(Frequency ~ InsLengthGroup, data = L1Ewing2011_S4, FUN = min)
FreqMax <- aggregate(Frequency ~ InsLengthGroup, data = L1Ewing2011_S4, FUN = max)
FreqSum <- aggregate(Frequency ~ InsLengthGroup, data = L1Ewing2011_S4, FUN = sum)
FreqVar  <- aggregate(Frequency ~ InsLengthGroup, data = L1Ewing2011_S4, FUN = var)
plot(InsLengthMid, FreqMean$Frequency, 
     xlab = "L1 insertion length [bp]",
     ylab = "Frequency")
AddErrorBars(MidX = InsLengthMid, MidY = FreqMean$Frequency, 
             ErrorRange = sqrt(FreqVar$Frequency / FreqSum$Frequency * 
                                 FreqMean$Frequency),
             TipWidth = 50)
CreateDisplayPdf("D:/L1polymORF/Figures/L1InsertionLengthVsFrequencyEwing2011.pdf")

# Plot quantiles of frequencies of small and large fragments against each other
FreqSmall <- L1Ewing2011_S4$Frequency[L1Ewing2011_S4$length <= 500] 
FreqMed <- L1Ewing2011_S4$Frequency[L1Ewing2011_S4$length >= 2000 &
                                        L1Ewing2011_S4$length <= 5000] 
FreqLarge <- L1Ewing2011_S4$Frequency[L1Ewing2011_S4$length >= 6000] 

# Create quantile plots
par(mfrow = c(1,2))
qqplot(FreqSmall, FreqMed, xlab = "Frequencies of small fragments",
       ylab = "Frequencies of medium fragments")
lines(c(0, 1), c(0, 1))
qqplot(FreqLarge, FreqMed, xlab = "Frequencies of large fragments",
       ylab = "Frequencies of medium fragments")
lines(c(0, 1), c(0, 1))
CreateDisplayPdf("D:/L1polymORF/Figures/L1InsertionLengthVsFrequencyQQPlots.pdf",
                 height = 4)
