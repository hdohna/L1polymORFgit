# The script below reads in read length data from two fastq files and plots 
# histograms to compare them

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load read length data
load("D:/L1polymORF/Data/PacBioReadLengths.RData")

# Create plot
HistPS     <- hist(RL1, plot = F, breaks = seq(0, 45000, 500))
HistNormal <- hist(RL2, plot = F, breaks = seq(0, 50000, 500))
plot(HistPS$mids - 50, HistPS$density, type = "s",
     col = "red", ylim = c(0, 0.0006), xlim = c(0, 15000),
     xlab = "Read length", ylab = "Frequency")
lines(HistNormal$mids + 50, HistNormal$density, type = "s", col = "blue")
legend(x = 10000, y = 5.5*10^-4, legend = c("Old polymerase", "New polymerase"), lty = c(1, 1),
       col = c("blue", "red"), bty = "n")
CreateDisplayPdf("D:/L1polymORF/Figures/ReadLDistns.pdf",
                 PdfProgramPath = '"C:\\Program Files (x86)\\Adobe\\Reader 11.0\\Reader\\AcroRd32"',
                 height = 7)

