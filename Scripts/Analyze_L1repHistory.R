# Load data about replication history (generated in script L1HS_repHistory.R)
load("D:/L1polymORF/Data/L1repHistoryResults.RData")

# Get table with L1 in reference genme from repeat masker
L1Table <- read.csv("D:/L1polymORF/Data/L1HS_repeat_table.csv", as.is = T)
L1HSSeqNames <- paste(L1Table$genoName, L1Table$genoStart, L1Table$genoEnd,
                      L1Table$strand, sep = "_")

# Match each column in CountMat to table entry
NameMatch <- match(colnames(CountMat), L1HSSeqNames)

# Get activity and total number of fragments
Act <- as.numeric(L1Table$Act_L1rp[NameMatch])
Freq <- as.numeric(L1Table$Allele_frequency.[NameMatch])
WeightAct <- Act*Freq
TotalFrags <- colSums(CountMat)

# Fit a linear model that relates activity to total number of fragments
DistVals <- as.numeric(rownames(CountMat))
MeanAge <- apply(CountMat, 2, function(x) (x %*% DistVals) / sum(x) )
Fit1 <- lm(Act ~ CountMat[1,] + TotalFrags + MeanAge)
summary(Fit1)
Fit2 <- lm(Act ~ CountMat[1,])
summary(Fit2)
Fit3 <- lm(CountMat[1,] ~ WeightAct)
summary(Fit3)

# Plot activity vs number of close distance fragments
plot(Act, CountMat[1,], xlab = "Activity", 
     ylab = "Number of close distance fragments")
CreateDisplayPdf("D:/L1polymORF/Figures/ActVsNrFragm.pdf")
text(CountMat[1,], Act,TotalFrags, pos = 3)
cor.test(CountMat[1,!is.na(Act)], Act[!is.na(Act)])
cor.test(MeanAge[!is.na(Act)], Act[!is.na(Act)])
plot(MeanAge, Act, xlab = "Mean fragment age", 
     ylab = "Activity", xlim = c(0, 0.1))
plot(WeightAct, CountMat[1,], xlab = "Weighted activity", 
     ylab = "Number of close distance fragments")

