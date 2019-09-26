# Load data with HWE values
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/HWE_L1_and Flanks.RData")
hist(HWE_L1$P,    breaks = seq(0, 1, 0.001))
hist(HWE_Left$P,  breaks = seq(0, 1, 0.001))
hist(HWE_Right$P, breaks = seq(0, 1, 0.001))

# Get observed and expected heterozygotes
ObsExpHet <- as.data.frame(t(sapply(1:nrow(HWE_L1), function(i){
  ObsHetChar <- strsplit(as.character(HWE_L1$OBS.HOM1.HET.HOM2.[i]), "/")[[1]][2]
  ExpHetChar <- strsplit(as.character(HWE_L1$E.HOM1.HET.HOM2.[i]), "/")[[1]][2]
  c(ObsHet = as.numeric(ObsHetChar), ExpHet = as.numeric(ExpHetChar))
})))

sum(HWE_L1$P < 0.05 & ObsExpHet$ObsHet > ObsExpHet$ExpHet)
sum(HWE_L1$P < 0.05)

# Make genomic ranges
HWEGR <- makeGRangesFromDataFrame(HWE_L1, seqnames.field = "CHR",
                                 start.field = "POS",
                                 end.field = "POS")

# Read repeat masker table for L1HS
L1Table <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1HS_repeat_table_Hg19.csv", as.is = T)
L1Table$ChrNr <- substr(L1Table$genoName, 4, nchar(L1Table$genoName))

# Create GRanges objects with L1 Seqences
L1GR <- makeGRangesFromDataFrame(L1Table, seqnames.field = "ChrNr",
                                 start.field = "genoStart",
                                 end.field = "genoEnd")

# Find overlap 
OL <- findOverlaps(HWEGR, L1GR)

HWEperL1 <- aggregate(HWE_L1$P[OL@to], by = list(OL@from), FUN = mean)
                      
plot(width(L1GR)[HWEperL1$Group.1], HWEperL1$x)
cor.test(width(L1GR)[HWEperL1$Group.1], HWEperL1$x)
cor.test(width(L1GR)[HWEperL1$Group.1], HWEperL1$x, method = "spearman")
