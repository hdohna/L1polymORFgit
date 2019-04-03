# The script below estimates selection coefficients of L1 from the 
# 1000 genome data

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

FreqVals <- seq(0, 1, 0.001)
SCoeffs  <- seq(-10^-4, 10^-4, 5*10^-5)
Cols     <- rainbow(length(SCoeffs))
Distn1 <- sapply(FreqVals, function(x) AlleleFreqDistn(y = x, s = 0, N = 10^5))
plot(FreqVals, Distn1, type = "l", ylim = c(0, 6),
     xlab = "Allele frequency",
     ylab = "Probability density")
for (i in 1:length(SCoeffs)){
  Distn <- sapply(FreqVals, function(x) AlleleFreqDistn(y = x, s = SCoeffs[i], N = 10^5))
  lines(FreqVals, Distn, col = Cols[i])
}
legend("topright", legend = SCoeffs, col = Cols, 
       lty = rep(1, length(SCoeffs)), y.intersp = 0.5)

CreateDisplayPdf('D:/L1polymORF/Figures/FreqDistPerSelectCoeff.pdf',
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7, width = 7)
