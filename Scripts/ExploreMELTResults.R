# The following script explores the results by Gardner et al. 2017 genome 
# reserach that were downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/
load('D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData')
load('D:/L1polymORF/Data/L1RefRanges_hg19.Rdata')

# Read in vcf file with variant calls
MEVarCall <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_call.vcf", 
                        as.is = T,
                        col.names = c("Chrom", "Pos", "ID", "Alt", "Type", "V6", 
                                      "V7", "Info"))
MEVarRegion <- read.table("D:/L1polymORF/Data/nstd144.GRCh37.variant_region.vcf", as.is = T)
table(MEVarRegion$V5)
MEVarRegion$V8[1]
# Read sample data
SampleInfo <- read.delim("D:/L1polymORF/Data/igsr_samples.tsv", as.is = T)
table(SampleInfo$Data.collections)
SampleInfo$Data.collections[1:5]

# Get index of LINE1 deletions and check whether an allele frequency is counted
idxDel <- grep("<DEL:ME:LINE1>", MEVarCall$Type)
MEVarCall$Info[idxDel[1:10]]
grep("AF=", MEVarCall$Info[idxDel])

idxIns <- grep("<INS:ME:LINE1>", MEVarCall$Type)
grep("AF=", MEVarCall$Info[idxIns])
length(idxIns)
idxIns[1]
MEVarCall$Info[idxIns[1]]

# Extract allele frequency from info column
GetAF <- function(x){
  xSplit <- strsplit(x, ";")[[1]]
  AFch   <- strsplit(xSplit[length(xSplit)], "=")[[1]][2]
  as.numeric(AFch)
}
MEVarCall$Freq <- sapply(MEVarCall$Info, GetAF)

# Create GRanges object for MEVarCall
MEVarCall$ChromName <- paste("chr", MEVarCall$Chrom, sep = "")
MEVar_GR <- makeGRangesFromDataFrame(df = MEVarCall,
                                     seqnames.field = "ChromName",
                                     start.field = "Pos",
                                     end.field = "Pos")

# Find overlap between the two GRanges and plot frequency estimates against 
# each other
OL_1000G_MEVar <- findOverlaps(L1_1000G_GR_hg19, MEVar_GR)
plot(L1_1000G_reduced$Frequency[OL_1000G_MEVar@from], 
     MEVarCall$Freq[OL_1000G_MEVar@to])
lines(c(0, 1), c(0, 1))

NewGR <- resize(MEVar_GR, 100, fix = "center")
sum(overlapsAny(L1GRanges,  NewGR))
OL_Ref_MEVar <- findOverlaps(L1GRanges,  MEVar_GR)
L1GRanges[OL_Ref_MEVar@from[1:10]]
MEVar_GR[OL_Ref_MEVar@to[1:10]]
