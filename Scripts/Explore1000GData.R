# The following script explores L1 insertion patterns

library(glmnet)

# Read in data
load("D:/L1polymORF/Data/GRanges_L1_1000Genomes.RData")
SampleInfo <- read.delim("D:/L1polymORF/Data/1000GenomeSampleInfo.txt")
SampleMatch <- match(SampleColumns, SampleInfo$sample)
SampleInfo <- SampleInfo[SampleMatch, ]

# Patterns of Nr L1 per individual
NrL1_perSample <- colSums(L1_1000G[,SampleColumns])
hist(NrL1_perSample, breaks = 0:300)

# Variation in L1 per individual by super population
NrL1BySuperPop <- lm(NrL1_perSample ~ SampleInfo$super_pop)
anova(NrL1BySuperPop)
boxplot(NrL1_perSample ~ SampleInfo$super_pop)

# Variation in L1 per individual by population
NrL1ByPop <- lm(NrL1_perSample ~ SampleInfo$pop)
anova(NrL1ByPop)
boxplot(NrL1_perSample ~ SampleInfo$pop, cex.axis = 0.5)

# Correlation between individual L1 and the total number of L1 
x <- 1
L1_TotCorr <- sapply(1:nrow(L1_1000G), function(x){
  cor(t(L1_1000G[x,SampleColumns]), NrL1_perSample)
})
hist(L1_TotCorr, breaks = seq(-0.2, 0.6, 0.001))
                     
# run Elastic net
ElasticFit <- cv.glmnet(ScoreMat, HISubset$logT, 
                        lower.limits = c(rep(0, NrT), rep(-Inf, ncol(ScoreMat) - NrT)))



