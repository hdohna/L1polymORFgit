# The script below estimates selection coefficients of L1 from the 
# 1000 genome data

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# Specify file paths
InputPath <- "D:/L1polymORF/Data/L1SelectionResults.RData"
InputPath <- "D:/L1polymORF/Data/L1SelectionResults_MELT.RData"

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Load previously generated objects
load(InputPath)

###################################################
#                                                 #
#   Plot density vs. selection coefficient        #
#                                                 #
###################################################

# Probability 
N1 <- length(L1GRanges)

# Create a vector of selection coefficients
SCoeffVect <- c(Promoter = ML_L1ExonIntron$par[1],
                Exon = sum(ML_L1ExonIntron$par[c(1, 2)]),
                Intron = sum(ML_L1ExonIntron$par[c(1, 3)]),
                Intergenic = ML_L1ExonIntron$par[1])
names(SCoeffVect) <- sapply(names(SCoeffVect), 
                            function(x) strsplit(x, "\\.")[[1]][1])

# Plot selection coefficient against 
if (!all(names(SCoeffVect) == colnames(InsPerbp))){
  stop("Selection coefficients and L1 densities are not in same order!")
}
if (!all(names(SCoeffVect) == names(MeanFreqs))){
  stop("Selection coefficients and L1 frequencies are not in same order!")
}

# Get sample size and create a range of s-values
SSize <- 2*2504
SVals <- seq(-0.0025, -0.00001, 0.00001)

par(oma = c(2, 1, 1, 3), mfcol = c(2, 2), mai = c(1, 1, 0.2, 0.2),
    cex.axis = 1, cex.lab = 1.5)
layout(rbind(1:2, c(3, 3)), widths = c(1, 1))
layout(rbind(c(1, 1, 2, 2), c(0, 3, 3, 0)), widths = c(1, 1))
# Plot LINE-1 frequency against number of LINE-1 per Mb
plot(InsPerbp[2,], MeanFreqs, xlab = "LINE-1s per Mb", 
     ylab = "Mean LINE-1 frequency", main = "A")
text(InsPerbp[2,] + 2*10^(-1)*c(1, 0, 0, -1.2), 
     MeanFreqs + c(3*10^-3, 3*10^-3, 3*10^-3, -3*10^-3),  
     names(SCoeffVect))

# Plot expected frequency versus observed mean frequency
ExpL1 <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = 2*2504))
plot(SCoeffVect, MeanFreqs, ylab = "Mean LINE-1 frequency", 
     xlab = "Selection coefficient", xlim = c(-0.0025, 0.0007), main = "B")
text(SCoeffVect + 2*c(0.0003, 0, -0.0003, -0.0003), 
     MeanFreqs + c(0, 3*10^-3, 0, 0), names(SCoeffVect))
lines(SVals, ExpL1)


# Plot probability for inclusion versus number of LINE-1 per Mb
ProbL1 <- sapply(SVals, function(x) ProbAlleleIncluded(x,N = 5*10^4, 
                                                       SampleSize = 2*2504))
plot(SCoeffVect, InsPerbp[2,], ylab = "LINE-1s per Mb", 
     xlab = "Selection coefficient", xlim = c(-0.0025, 0), ylim = c(0, 3),
     main = "C")
text(SCoeffVect, InsPerbp[2,] + 2*10^(-1), names(SCoeffVect))
par(new = TRUE)
plot(SVals, ProbL1, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",
     ylim = c(0, 1))
axis(side = 4)
mtext("Inclusion probability", 4, line = 3)
#mtext("Selection coefficient", 1, line = 3)

CreateDisplayPdf('D:/L1polymORF/Figures/SelectionPerRegion_MELT.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Plot frequency vs. insertion length           #
#                                                 #
###################################################

# Create a vector of L1 start classes
L1_1000G$L1StartNumClass <- cut(L1_1000G$L1StartNum, breaks = 
                                  seq(0, 6100, 500))
L1_1000G$Frequency

# Get mean L1 frequency per start
L1StartAggregated <- aggregate(L1_1000G[,c("L1StartNum", "Frequency")], 
                               by = list(L1_1000G$L1StartNumClass), FUN = mean)

# Get sample size and create a range of s-values
SSize <- 2*2504
StartVals  <- seq(0, 6000, 100)
Full       <- StartVals == 0
SVals <- ML_L1startL1full$par[1] + ML_L1startL1full$par[2]*StartVals +
  ML_L1startL1full$par[3]*Full

# Plot expected frequency versus observed mean frequency
ExpL1Start <- sapply(SVals, function(x) ExpAlleleFreq(x, N = 10^4, SampleSize = 2*2504))
par( mfrow = c(1, 1))
plot(L1StartAggregated$L1StartNum, 
     L1StartAggregated$Frequency * SSize, xlab = "5' start of LINE-1",
     ylab = "Mean LINE-1 frequency")
lines(StartVals, ExpL1Start)
mtext("Selection coefficient", 1, line = 3)
CreateDisplayPdf('D:/L1polymORF/Figures/FreqVsL1Start.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)

###################################################
#                                                 #
#   Fit effect of strandedness on selection       #
#                                                 #
###################################################

# Create a matrix of predictor variables (L1 start and boolean variable for)
PredictMatWithinGene <- L1_1000G[L1_1000G$blnOLGene, 
                                 c( "blnOLGeneSameStrand", "blnOLGene", "blnOLGene")]

# Estimate maximum likelihood for a single selection coefficient
ML_1Par_gene <- constrOptim(theta = c(a = -0.0004),
                            f = function(x) -AlleleFreqLogLik_4Par(
                              Freqs = (L1_1000G$Frequency[L1_1000G$blnOLGene] * 2*2504),
                              Counts = rep(1, sum(L1_1000G$blnOLGene)),
                              Predict = PredictMatWithinGene,
                              a = x[1], b = 0, c = 0, d = 0, N = 10^4,
                              SampleSize = 2*2504),
                            grad = NULL,
                            ui = rbind(1,-1),
                            ci = c(a = -0.01, a = -0.01),
                            method = "Nelder-Mead")

# Get maximum likelihood estimate for effect of exonic L1 on selection
cat("Estimate effect of same strand overlap on selections ...")
ML_L1SameStrand <-  constrOptim(theta = c(a = ML_1Par_gene$par, b = 0),
                          f = function(x) -AlleleFreqLogLik_4Par(
                            Freqs = (L1_1000G$Frequency[L1_1000G$blnOLGene] * 2*2504),
                            Counts = rep(1, sum(L1_1000G$blnOLGene)),
                            Predict = PredictMatWithinGene,
                            a = x[1], b = x[2], c = 0, d = 0, N = 10^4,
                            SampleSize = 2*2504),
                          grad = NULL,
                          ui = rbind(c(1, 0),  c(0, 1),   
                                     c(-1, 0), c(0, -1)),
                          ci = c(a = -0.01, b = -10^(-2), 
                                 a = -0.01, b = -10^(-2)),
                          method = "Nelder-Mead")
cat("done!\n")

# Get columns of AIC and parameter values
Cols2Append <- t(sapply(list(ML_1Par_gene, ML_L1SameStrand), 
                        function(x){
                          c(NrParameters = GetNPar(x), AIC = GetAIC(x), 
                            Pars = GetParVals(x))
                        }))
# Combine AIC values into one vector
AICTabWithinGene <- cbind(data.frame(
  Predictor = c("none", "SameStrand"),
  stringsAsFactors = F),
  Cols2Append)


###################################################
#                                                 #
#   Fit effect of singleton coef. on selection    #
#                                                 #
###################################################

# Create a matrix of predictor variables 
PredictMat <- L1SingletonCoeffs[, c("coef", "L1Start", "L1Start")]
blnNA <- sapply(1:nrow(L1SingletonCoeffs), function(x) any(is.na(PredictMat[x,])))

# Determine maximum likelihood with one parameter (selection coefficient)
cat("Maximizing likelihood for one parameter (selection coefficient) ...")
ML_1Par_coef <- optim(par = c(a = ML_1Par$par),
                      fn = function(x) -AlleleFreqLogLik_4Par(
                        Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA], 
                        Counts = rep(1, sum(!blnNA)), 
                        Predict = PredictMat[!blnNA,], 
                        a = x, b = 0, c = 0, d = 0, N = 10^4, 
                        SampleSize = 2*2504),
                      lower = -0.01, upper = 0.01,
                      method = "L-BFGS-B")
ML_2Pars_L1coef <- constrOptim(
  theta = c(a = ML_1Par_coef$par, b = 0),
  f = function(x) -AlleleFreqLogLik_4Par(
    Freqs = (L1SingletonCoeffs$Freq * 2*2504)[!blnNA], 
    Counts = rep(1, sum(!blnNA)), 
    Predict = PredictMat[!blnNA,], 
    a = x[1], b = x[2], c = 0, d = 0, N = 10^4, 
    SampleSize = 2*2504),
  grad = NULL,
  ui = rbind(c(1, 0),  c(0, 1),  
             c(-1, 0), c(0, -1)),
  ci = c(a = -0.01, b = -2*10^(-3), 
         a = -0.01, b = -2*10^(-3)),
  method = "Nelder-Mead")
cat("done!\n")

# Save everything
save.image(SelectResultOutPath)
