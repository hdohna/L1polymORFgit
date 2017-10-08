# The following script analyzes the correlation between L1 frequency and 
# insertion length in 1000 genome data

# Load packages
library(ggplot2)
source("D:/L1polymORFgit/Scripts/_Start_L1polymORF.R")

# Load L1 catalog GenomicRanges
load("D:/L1polymORF/Data/L1CatalogGRanges.RData")
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')

# Median frequency of fragments and full-length L1
median(L1_1000G_reduced$Frequency[L1_1000G_reduced$InsLength <= 5900], na.rm = T)
median(L1_1000G_reduced$Frequency[L1_1000G_reduced$InsLength > 6000], na.rm = T)

# Correlation between frequency and insertion length among fragments
cor.test(L1_1000G_reduced$Frequency[L1_1000G_reduced$InsLength <= 5900], 
         L1_1000G_reduced$InsLength[L1_1000G_reduced$InsLength <= 5900], method = "spearman")
plot(L1_1000G_reduced$Frequency, L1_1000G_reduced$InsLength)
cor.test(L1_1000G_reduced$Frequency[L1_1000G_reduced$InsLength > 6000], 
         L1_1000G_reduced$InsLength[L1_1000G_reduced$InsLength > 6000], method = "spearman")
ggplot(L1_1000G_reduced, aes(x = InsLength, y = Frequency)) + 
  geom_point(alpha = 0.1, size = 3) +
  geom_smooth()
  
# Group insertion length in units of 500 bp
InsLengthGroups <- seq(0, 6500, 500)
InsLengthMid <- InsLengthGroups[-length(InsLengthGroups)] + 250
L1_1000G_reduced$InsLengthGroup <- cut(L1_1000G_reduced$InsLength,
                                       breaks = InsLengthGroups)
FreqMean <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = mean)
FreqSum <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = sum)
FreqVar  <- aggregate(Frequency ~ InsLengthGroup, data = L1_1000G_reduced, FUN = var)
plot(InsLengthMid, FreqMean$Frequency, ylim = c(0, 0.05), 
     xlab = "L1 insertion length [bp]",
     ylab = "Frequency")
AddErrorBars(MidX = InsLengthMid, MidY = FreqMean$Frequency, 
             ErrorRange = sqrt(FreqVar$Frequency / FreqSum$Frequency * 
                                 FreqMean$Frequency),
             TipWidth = 50)
CreateDisplayPdf("D:/L1polymORF/Figures/L1InsertionLengthVsFrequency1000Genomes.pdf")
