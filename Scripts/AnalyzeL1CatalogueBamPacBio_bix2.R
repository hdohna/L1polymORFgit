# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and determines which L1 is present

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
# source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)

# Specify parameters
BorderWidth <- 100

# Specify file paths
BamFile <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue.dedup.unique.sorted.bam'
#BamFile <- '/share/diskarray3/hzudohna/NA12878PacBio_aln2Catalogue2016-05-07filteredByLength.bam'
CatalogueFile <- '/home/hzudohna/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv'
CatalogSeqFile <- "/home/hzudohna/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fas"
AlignListFile <- "/home/hzudohna/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016_L1Locations.RData"
OutFile      <- "/home/hzudohna/L1polymORF/Data/L1CatalogPacBioSummary.RData"
############################
#                          #
#        Read Data         #
#                          #
############################

# Load file with  
load(file = AlignListFile)
if(!all(colnames(L1StartEnd) == names(L1withFlank))){
  stop("Names L1StartEnd and L1withFlank don't match!")
}
# Get accession numbers of sequences that contain L1 
blnSeqWithL1     <- L1StartEnd["Width", ] >= 5000
all(blnSeqWithL1)

# Read catalog table
L1Catalog <- read.csv(CatalogueFile, as.is = T)
AccMatch <- match(colnames(L1StartEnd), L1Catalog$Accession)
L1CatalogMatched <- L1Catalog[AccMatch, ] 

# Get genomic ranges of catalogue
GRL1Catalogue <- GRanges(seqnames = colnames(L1StartEnd), 
                         ranges = IRanges(start = L1StartEnd["Start", ],
                                          end = L1StartEnd["End", ]),
                         strand = L1CatalogMatched$Strand)

# Get reads per catalogue entry
cat("********    Extracting reads   ******************\n\n")
ReadList <- lapply(1:length(GRL1Catalogue), function(i){
  Chrom       <- colnames(L1StartEnd)[i]
  ChromLength <- length(L1withFlank[[i]])
  R1 <- GRL1Catalogue[i]
  cat("Extracting reads of L1", Chrom, "\n")
  Reads <- extractReads(bam.file = BamFile, region = R1)
})
names(ReadList) <- names(L1withFlank)

# Get number of reads per catalogue entry
NrReads <- sapply(ReadList, length, USE.NAMES = T)
pdf("/home/hzudohna/L1polymORF/Figures/Na12878Pacbio_CatalogueReadHist.pdf")
hist(NrReads, xlab = "Number of aligned reads", ylab = "Number of LINE1")
dev.off()

################################
#                              #
#     Analyze coverage         #
#                              #
################################

cat("\n********    Summarizing coverage   ******************\n\n")

# Create a list of coverage vectors
CoverMat<- t(sapply(ReadList, function(x){
  Cover <- coverage(x, width = 16000)
  as.vector(Cover[[1]])
}))

# Invert coverage row vectors for L1 on the negative strand
idxNegStrand <- which(L1CatalogMatched$Strand == "-")
for (i in idxNegStrand){
  CoverMat[i, ] <- CoverMat[i, 16000:1]
}

# Get meadians and 95% quantiles 
QuantileMat <- apply(CoverMat, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
pdf("/home/hzudohna/L1polymORF/Figures/Na12878Pacbio_CatalogueAverageCover.pdf")
plot(QuantileMat[2,], type = "n", ylim = c(0, 5000), 
     ylab = 'Coverage', xlab = "Genomic position")
idxFw <- 1:ncol(CoverMat)
idxRv <- ncol(CoverMat):1
polygon(c(idxFw, idxRv), c(QuantileMat[1, idxFw], QuantileMat[3, idxRv]),
        col = "grey", border = NA)
lines(QuantileMat[2,])
segments(c(5000, 11000), c(0,0), c(5000, 11000), c(10^5, 10^5), col = "red",
         lty = 2)
dev.off()

################################
#                              #
#     Filter reads             #
#                              #
################################

cat("\n********    Filtering reads   ******************\n\n")

# Get reads per catalogue entry
ReadListFiltered <- lapply(1:length(ReadList), function(i){
  FRange <- GRanges(seqnames = colnames(L1StartEnd)[i], 
                    IRanges(start = 1, end = 16000))
  L1Borders <- GRanges(seqnames = colnames(L1StartEnd)[i], 
                       IRanges(start = L1StartEnd[1:2,i] - BorderWidth, 
                               end = L1StartEnd[1:2,i] + BorderWidth))
  subsetByOverlaps(Reads, L1Borders, minoverlap = 2 * BorderWidth)
})
names(ReadList) <- names(L1withFlank)

# Get number of reads per catalogue entry
NrReadsFiltered <- sapply(ReadListFiltered, length, USE.NAMES = T)
pdf("/home/hzudohna/L1polymORF/Figures/Na12878Pacbio_CatalogueReadHist_Filtered.pdf")
hist(NrReadsFiltered, xlab = "Number of aligned reads", ylab = "Number of LINE1")
dev.off()

##################################################
#                                                #
#     Analyze coverage of filtered reads         #
#                                                #
##################################################

cat("\n********   Summarizing coverage   ******************\n\n")

# Create a list of coverage vectors
CoverMatFiltered<- t(sapply(ReadListFiltered, function(x){
  Cover <- coverage(x, width = 16000)
  as.vector(Cover[[1]])
}))

# Invert coverage row vectors for L1 on the negative strand
idxNegStrand <- which(L1CatalogMatched$Strand == "-")
for (i in idxNegStrand){
  CoverMatFiltered[i, ] <- CoverMatFiltered[i, 16000:1]
}

# Get meadians and 95% quantiles 
QuantileMatFiltered <- apply(CoverMatFiltered, 2, FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
pdf("/home/hzudohna/L1polymORF/Figures/Na12878Pacbio_CatalogueAverageCoverFiltered.pdf")
plot(colMeans(CoverMatFiltered), type = "l", 
     xlab = "Genomic position", ylab = "Coverage")
segments(c(5000, 11000), c(0,0), c(5000, 11000), c(10^5, 10^5), col = "red",
         lty = 2)
dev.off()
pdf("/home/hzudohna/L1polymORF/Figures/Na12878Pacbio_CatalogueAverageCoverFiltered_lowRead.pdf")
plot(colMeans(CoverMatFiltered[NrReadsFiltered < 10, ]), type = "l", 
     xlab = "Genomic position", ylab = "Coverage")
dev.off()


# # Table combination of large number of reads and presenc in reference
# L1CatalogueMatched <- L1Catalogue[match(names(NrReads), L1Catalogue$Accession), ]
# blnReference <- abs(L1CatalogueMatched$end_HG38 - L1CatalogueMatched$start_HG38) > 5000
# blnManyReads <- NrReads > 1000
# table(blnReference, blnManyReads)
# 
# # Save accession numbers of L1s with reads
# write.csv(L1CatalogueMatched[blnManyReads, ],
#           "/home/hzudohna/L1polymORF/Data/L1CatalogueWithReads_PacBio.csv")
# 

cat("\n********   Saving results   ******************\n\n")
save(list = c("ReadList", "NrReads", "CoverMat", "QuantileMat",
              "ReadListFiltered", "NrReadsFiltered", "CoverMatFiltered", 
              "QuantileMatFiltered"), OutFile)
cat("Results saved to", OutFile)


