# Explore the 1000 genome vcf file found at 
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/paper_data_sets/companion_papers/mapping_structural_variation/
# from Mills et al. 2011 Nature: http://www.nature.com/nature/journal/v470/n7332/full/nature09708.html

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
names(BSgenome.Hsapiens.UCSC.hg38)

# Read in table with known L1 
L1Catalogue <- read.csv("D:/L1polymORF/Data/L1Catalogue_Fri_Apr_01_18-28-08_2016.csv", 
                        as.is = T)

# Read in table with L1 uhn et al. 2014
L1Kuhn2014 <- read.csv("D:/L1polymORF/Data/Kuhn_2014_Table_S2.csv", 
                        as.is = T)

# Read vcf file from 1000 Genome data
MEI1000Gvcf <- ReadVCF("D:/L1polymORF/Data/union.2010_06.MobileElementInsertions.genotypes.vcf")

# Subset to obtain only L1s
L1Ins1000G <- MEI1000Gvcf[MEI1000Gvcf$ALT == "<INS:ME:L1>", ]

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

# Determine sample columns and get genotypes
SampleColumns <- grep("NA", colnames(L1Ins1000G), value = T)
GenoDF <- sapply(SampleColumns, function(x)GetGenotype(L1Ins1000G[,x]))
colnames(GenoDF) <- paste("Geno", SampleColumns, sep = "_")
L1Ins1000Ggeno <- cbind(L1Ins1000G, GenoDF)
L1Ins1000G$InsLength > 6000

# Explore genotype data
hist(rowSums(GenoDF[L1Ins1000G$InsLength > 6000, ], na.rm = T)/ 2/rowSums(GenoWithData[L1Ins1000G$InsLength > 6000, ]))
hist(colSums(GenoDF, na.rm = T))
hist(colSums(GenoDF > 0, na.rm = T) / colSums(GenoWithData),
     xlab = "Proportion of L1 insertions present")
hist(rowSums(GenoDF, na.rm = T) / rowSums(GenoWithData)) 

# Lift 1000Genome coordinates over to Hg38 
GRL1Ins1000G <- GRanges(seqnames = paste("chr", L1Ins1000Ggeno$X.CHROM, sep = ""),
                        ranges = IRanges(start = L1Ins1000Ggeno$POS,  
                                         end = L1Ins1000Ggeno$POS))
GRL1Ins1000G_hg38    <- liftOver(GRL1Ins1000G, 
                           chain = import.chain("D:/L1polymORF/Data/hg18ToHg38.over.chain"))
NrMapped_hg38        <- sapply(GRL1Ins1000G_hg38, length)
idxUniqueMapped_hg38 <- which(NrMapped_hg38 == 1) 
GRL1Ins1000G_hg38Mapped <- unlist(GRL1Ins1000G_hg38[idxUniqueMapped_hg38])

# Lift Kuhn coordinates over to Hg38 
GRL1Kuhn2014 <- GRanges(seqnames = L1Kuhn2014$chromosome,
                        ranges = IRanges(start = L1Kuhn2014$start,  
                                         end = L1Kuhn2014$end),
                        strand = L1Kuhn2014$strand)
GRL1Kuhn2014_hg38    <- liftOver(GRL1Kuhn2014, 
                                 chain = import.chain("D:/L1polymORF/Data/hg19ToHg38.over.chain"))
NrMapped_hg38        <- sapply(GRL1Kuhn2014_hg38, length)
idxUniqueMapped_hg38 <- which(NrMapped_hg38 == 1) 
GRL1Kuhn2014_hg38Mapped <- unlist(GRL1Kuhn2014_hg38[idxUniqueMapped_hg38])


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
hist(Dists, breaks = seq(0, 1.5*10^8, 10^5))
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

