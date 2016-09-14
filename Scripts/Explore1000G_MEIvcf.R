# Explore the 1000 genome vcf file found at 
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/paper_data_sets/companion_papers/mapping_structural_variation/
# from Mills et al. 2011 Nature: http://www.nature.com/nature/journal/v470/n7332/full/nature09708.html

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
names(BSgenome.Hsapiens.UCSC.hg38)

# Specify file paths
GROutputPath <- 'D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData'

############################
#                          #
#        Read Data         #
#                          #
############################

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalogue_Fri_Apr_01_18-28-08_2016.csv", 
                        as.is = T)

# Read in table with L1 Kuhn et al. 2014
L1Kuhn2014 <- read.csv("D:/L1polymORF/Data/Kuhn_2014_Table_S2.csv", as.is = T)
Kuhn2014Samples <- read.csv("D:/L1polymORF/Data/Kuhn et al 2014 PNAS Table S1.csv", 
                       as.is = T)
Kuhn2014Genotypes <- read.csv("D:/L1polymORF/Data/Kuhn et al 2014 PNAS Table S4.csv", 
                                   as.is = T)
Kuhn2014Primers <- read.csv("D:/L1polymORF/Data/Kuhn et al 2014 PNAS Table S5.csv", 
                              as.is = T)
Kuhn2014Genotypes1000G <- read.csv("D:/L1polymORF/Data/Kuhn et al 2014 PNAS Table S6.csv", 
                            as.is = T)

# Read vcf file from 1000 Genome data
MEI1000Gvcf <- ReadVCF("D:/L1polymORF/Data/union.2010_06.MobileElementInsertions.genotypes.vcf")
MEI1000GLines <- readLines("D:/L1polymORF/Data/union.2010_06.MobileElementInsertions.genotypes.vcf")
MEI1000GLines[1:32]

# Look up which names in Kuhn at al 2014 are in 1000 genome dataset
colnames(Kuhn2014Genotypes1000G)[colnames(Kuhn2014Genotypes1000G) %in% MEI1000Gvcf$ID]
colnames(Kuhn2014Genotypes1000G)[!colnames(Kuhn2014Genotypes1000G) %in% MEI1000Gvcf$ID]

# Subset to obtain only L1s
L1Ins1000G <- MEI1000Gvcf[MEI1000Gvcf$ALT == "<INS:ME:L1>", ]

##############################################
#                                            #
#   Get insertion length and genotype        #
#                                            #
##############################################

# Get insertion length
L1Ins1000G$InsLength <- sapply(L1Ins1000G$INFO, function(x){
  Xsplit <- strsplit(x, ";")[[1]]
  SVLENchar <- grep("SVLEN", Xsplit, value = T)
  if (length(SVLENchar) > 0){
    as.numeric(substr(SVLENchar, 7, nchar(SVLENchar)))
  } else {
    NA
  }
})

# Explore insertion length
hist(L1Ins1000G$InsLength)
sum(L1Ins1000G$InsLength > 6000, na.rm = T)

# function to get the genotype
GetGenotype <- function(SampleColumn){
  sapply(SampleColumn, function(x){
    Split1 <- strsplit(x, ":")[[1]][1]
    Split2 <- strsplit(Split1, "/")[[1]]
    sum(as.numeric(Split2))
  })
}
Geno <- GetGenotype(L1Ins1000G$NA06986)
table(Geno)

# Determine sample columns and get genotypes
SampleColumns <- grep("NA", colnames(L1Ins1000G), value = T)
GenoDF <- sapply(SampleColumns, function(x)GetGenotype(L1Ins1000G[,x]))
colnames(GenoDF) <- paste("Geno", SampleColumns, sep = "_")
L1Ins1000Ggeno <- cbind(L1Ins1000G, GenoDF)
L1Ins1000G$InsLength > 6000
GenoWithData <- !is.na(GenoDF)
sum(L1Ins1000G$InsLength > 6000, na.rm = T)

# Explore genotype data
L1Ins1000G$Freq <- rowSums(GenoDF, na.rm = T)
hist(rowSums(GenoDF[L1Ins1000G$InsLength > 6000, ], na.rm = T)/ 2 /
       rowSums(GenoWithData[L1Ins1000G$InsLength > 6000, ]))
hist(colSums(GenoDF, na.rm = T))
hist(colSums(GenoDF > 0, na.rm = T) / colSums(GenoWithData),
     xlab = "Proportion of L1 insertions present")
hist(rowSums(GenoDF, na.rm = T) / rowSums(GenoWithData)) 

# Create info dataset without the genotype columns
L1Ins1000G_Info <- L1Ins1000G[,!colnames(L1Ins1000G) %in% SampleColumns]

##############################################
#                                            #
#   Compare genomic coordinates              #
#                                            #
##############################################

# Lift 1000Genome coordinates over to Hg38 
GRL1Ins1000G <- GRanges(seqnames = paste("chr", L1Ins1000Ggeno$X.CHROM, sep = ""),
                        ranges = IRanges(start = L1Ins1000Ggeno$POS,  
                                         end = L1Ins1000Ggeno$POS))
GRL1Ins1000G_hg38    <- liftOver(GRL1Ins1000G, 
                           chain = import.chain("D:/L1polymORF/Data/hg18ToHg38.over.chain"))
length(GRL1Ins1000G_hg38)
NrMapped_hg38        <- sapply(GRL1Ins1000G_hg38, length)
length(NrMapped_hg38)
idxUniqueMapped_hg38 <- which(NrMapped_hg38 == 1) 
GRL1Ins1000G_hg38Mapped <- unlist(GRL1Ins1000G_hg38[idxUniqueMapped_hg38])
length(idxUniqueMapped_hg38)

# Lift Kuhn coordinates over to Hg38 
Kuhn2014PrimersNotNA <- Kuhn2014Primers[!is.na(Kuhn2014Primers$left.coordinate),]
GRKuhn2014Primers <- GRanges(seqnames = Kuhn2014PrimersNotNA$chromosome,
                        ranges = IRanges(start = Kuhn2014PrimersNotNA$left.coordinate,  
                                         end = Kuhn2014PrimersNotNA$left.coordinate),
                        strand = Kuhn2014PrimersNotNA$strand)
GRL1Kuhn2014 <- GRanges(seqnames = L1Kuhn2014$chromosome,
                        ranges = IRanges(start = L1Kuhn2014$start,  
                                         end = L1Kuhn2014$end),
                        strand = L1Kuhn2014$strand)
max(width(GRL1Kuhn2014))
overlapsAny(GRL1Kuhn2014, GRKuhn2014Primers)
GRL1Kuhn2014_hg38    <- liftOver(GRL1Kuhn2014, 
                                 chain = import.chain("D:/L1polymORF/Data/hg19ToHg38.over.chain"))
NrMappedKuhn_hg38        <- sapply(GRL1Kuhn2014_hg38, length)
idxUniqueMappedKuhn_hg38 <- which(NrMappedKuhn_hg38 == 1) 
GRL1Kuhn2014_hg38Mapped <- unlist(GRL1Kuhn2014_hg38[idxUniqueMappedKuhn_hg38])


# Create genomic ranges from catalogue
L1CatalogueMapped <- L1Catalogue[!is.na(L1Catalogue$start_HG38),]
GRCatalogue <- GRanges(seqnames = L1CatalogueMapped$Chromosome,
                        ranges = IRanges(start = L1CatalogueMapped$start_HG38,  
                                         end = L1CatalogueMapped$end_HG38),
                      strand = L1CatalogueMapped$Strand)

# Get for each mapped L1 insertion from Kuhn 2014 the distance to the 
# nearest L1 in the catalogue 
GRDistKuhn2014 <- distanceToNearest(GRL1Kuhn2014_hg38Mapped, GRCatalogue)
DistsKuhn2014  <- GRDistKuhn2014@elementMetadata@listData$distance
min(DistsKuhn2014)
DistsKuhn2014[order(DistsKuhn2014)[1:10]]
GRDist_Kuhn2014_1000G <- distanceToNearest(GRL1Kuhn2014_hg38Mapped, GRL1Ins1000G_hg38Mapped)
Dists_Kuhn2014_1000G  <- GRDist_Kuhn2014_1000G@elementMetadata@listData$distance
Dists_Kuhn2014_1000G[order(Dists_Kuhn2014_1000G)[1:1000]]

# Get for each mapped L1 insertion from 1000 genomes the distance to the 
# nearest L1 in the catalogue 
GRDist <- distanceToNearest(GRL1Ins1000G_hg38Mapped, GRCatalogue)
Dists <- GRDist@elementMetadata@listData$distance
min(Dists)
Dists[order(Dists)[1:10]]
hist(Dists, breaks = seq(0, 1.5*10^8, 10^5))


##############################################
#                                            #
#   Save genomic coordinates                 #
#                                            #
##############################################

save(list = c("GRL1Ins1000G_hg38Mapped", "idxUniqueMapped_hg38", "GRL1Ins1000G", 
              "L1Ins1000G_Info"), 
     file = GROutputPath)


##############################################
#                                            #
#   Match genotypes and L1 from Kuhn 2014    #
#                                            #
##############################################

# Create merged dataset
colnames(Kuhn2014Genotypes)[colnames(Kuhn2014Genotypes) == "X"] <- "L1.name"
KuhnL1GenoMerged <- merge(L1Kuhn2014, Kuhn2014Genotypes)
KuhnL1GenoMerged$Length <- KuhnL1GenoMerged$end - KuhnL1GenoMerged$start

# Look at L1 that are unique to one dataset
setdiff(Kuhn2014Genotypes$L1.name, L1Kuhn2014$L1.name)
setdiff(L1Kuhn2014$L1.name, Kuhn2014Genotypes$L1.name)

##############################################
#                                            #
#   Compare expected and observed freqs      #
#                                            #
##############################################

# Get L1 insertion frequencies
Freq <- rowSums(GenoDF, na.rm = T) / 2 / rowSums(GenoWithData, na.rm = T)

# Indicator for full-length insertion
blnFullLength <- L1Ins1000G$InsLength > 6000
hist(Freq[blnFullLength])
hist(Freq[!blnFullLength])
mean(Freq[blnFullLength], na.rm = T)
mean(Freq[!blnFullLength], na.rm = T)

# Get genotype sums
GenoMeans <- colSums(GenoDF, na.rm = T) / 2 / colSums(GenoWithData, na.rm = T)
GenoMeansFullL <- colSums(GenoDF[blnFullLength, ], na.rm = T) / 2 / 
  colSums(GenoWithData[blnFullLength, ], na.rm = T)
GenoMeansFragm <- colSums(GenoDF[!blnFullLength, ], na.rm = T) / 2 / 
  colSums(GenoWithData[!blnFullLength, ], na.rm = T)
hist(GenoMeans)
hist(GenoMeansFullL)
hist(GenoMeansFragm)
length(Freq)
dim(GenoDF)
rbinom(sum(!is.na(Freq)), 2, Freq[!is.na(Freq)])
x   <- 1
ind <- GenoDF[,x]
sum(rbinom(sum(L1Sampled), 2, Freq[L1Sampled]))

# Function to create qq-plot for a subset of L1 insertions
L1InsQQPlot <- function(idx, NrSamples = 1000){
  GenoWithDataFL <- GenoWithData[idx, ]
  FreqFL <- Freq[idx]
  SampleMat <- sapply(1:NrSamples, function(y){
    sapply(1:ncol(GenoDF), function(x) {
      L1Sampled <- GenoWithDataFL[,x] 
      sum(rbinom(sum(L1Sampled), 2, FreqFL[L1Sampled]))
    })
  })
  qqplot(as.vector(SampleMat), colSums(GenoDF[idx,], na.rm = T),
         xlab = "Predicted", ylab = "Observed")
  lines(c(0, 100), c(0, 100))
}
L1InsQQPlot(which(blnFullLength))
L1InsQQPlot(which(!blnFullLength))
