# The following script compares singleton density coefficients that were 
# calculated using different methods (phased and unphased)

# Load the first dataset
load('D:/L1polymORF/Data/SingletonAnalysis_unphased.RData')

# Rename the singleton coefficient dataframe
L1SingletonCoeffs_unphased <- L1SingletonCoeffs

# Load the second dataset
load('D:/L1polymORF/Data/SingletonAnalysis_phasedNoDir.RData')

# Match the two singleton coefficient datasets
ChrPos_unphased <- paste(L1SingletonCoeffs_unphased$Chrom, 
                         L1SingletonCoeffs_unphased$Pos, sep = "_")
ChrPos_phased <- paste(L1SingletonCoeffs$Chrom, 
                         L1SingletonCoeffs$Pos, sep = "_")
ChrPosMatch <- match(ChrPos_unphased, ChrPos_phased)
L1SingletonCoeffs <- L1SingletonCoeffs[ChrPosMatch,]

# Plot coefficients against each other
plot(L1SingletonCoeffs_unphased$coef, L1SingletonCoeffs$coef)
lines(c(-10, 10), c(-10, 10), col = "red")
cor(L1SingletonCoeffs_unphased$coef, L1SingletonCoeffs$coef)

# Plot p-values against each other
plot(log10(L1SingletonCoeffs_unphased$Pr...z..), 
     log10(L1SingletonCoeffs$Pr...z..))
lines(c(-10, 10), c(-10, 10), col = "red")

# Load iHH dataset
load("D:/L1polymORF/Data/iHHScores.RData")
iHHScores <- iHHScores[!is.na(iHHScores$iHH),]
hist(iHHScores$iHH, breaks = seq(-4, 6, 0.1))

# Match iHH data to 1000 genome data
IDmatch <- match(iHHScores$L1ID, L1_1000G$ID)
iHHScores$Chrom <- L1_1000G$CHROM[IDmatch]
iHHScores$POS   <- L1_1000G$POS[IDmatch]
iHHScores$Freq  <- L1_1000G_reduced$Frequency[IDmatch]

# Compare iHH score and frequency
plot(iHHScores$Freq, iHHScores$iHH)
cor.test(iHHScores$Freq, iHHScores$iHH)

# Match the two singleton coefficient datasets to the iHH data
ChrPos_iHH      <- paste(iHHScores$Chrom, iHHScores$POS, sep = "_")
ChrPosMatch_iHH <- match(ChrPos_iHH, ChrPos_unphased)
iHHScores$SinglCoeff_unphased <- L1SingletonCoeffs_unphased$coef[ChrPosMatch_iHH]
iHHScores$SinglCoeff_phased <- L1SingletonCoeffs$coef[ChrPosMatch_iHH]
plot(iHHScores$SinglCoeff_unphased, iHHScores$iHH)
plot(iHHScores$SinglCoeff_phased, iHHScores$iHH)
cor.test(iHHScores$SinglCoeff_unphased, iHHScores$iHH)
cor.test(iHHScores$SinglCoeff_phased, iHHScores$iHH)
