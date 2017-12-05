# The following script reads in data by Ewing et al (2011, Genome Research) and
# and combines them

# Source start script
source('D:/L1polymORFgit/Scripts/_Start_L1polymORF.R')

# Load packages
library(GenomicRanges)
library(rtracklayer)

# Read in both supplementarty tables by Ewing et al. 2011
L1Ewing2011_S2 <- read.csv("D:/L1polymORF/Data/Ewing2011Table_S2.csv")
L1Ewing2011_S4 <- read.csv("D:/L1polymORF/Data/Ewing2011Table_S4.csv")
L1Ewing2011_S4$chromosome <- paste("chr", L1Ewing2011_S4$Chr, sep = "")
L1Ewing2011_S2$chromosome <- paste("chr", L1Ewing2011$Chr, sep = "")
L1Ewing2011_S4_GR <- makeGRangesFromDataFrame(L1Ewing2011_S4, ignore.strand = T,
                                              seqnames.field = "chromosome",
                                              start.field = "Window.start",
                                              end.field = "Window.stop")
L1Ewing2011_S2_GR <- makeGRangesFromDataFrame(L1Ewing2011_S2, ignore.strand = T,
                                           seqnames.field = "chromosome",
                                           start.field = "Avg..Loc.",
                                           end.field = "Avg..Loc.")
findOverlaps(L1Ewing2011_S2_GR, L1Ewing2011_S4_GR)
