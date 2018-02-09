# The following scripts reads in genotype data from Seleme ets al. 2006 PNAS
# and analyzes whether genotypes with high L1 activity occur less frequently
# than predicted based on allele frequencies

# Read in genotype data from Seleme et al. 2006
GenoFreqs <- read.csv("D:/L1polymORF/Data/Seleme 2006 genotypes.csv", 
                      as.is = T)
GenoFreqs$isFemale <- 1:nrow(GenoFreqs) %in% 
  grep("female", GenoFreqs$Activity.categories)

# Create vectors of numeric activity sum values
ActivityNumSum <- as.numeric(as.character(GenoFreqs$ActivitySum))

# Retain only samples with numeric activity sum values 
GenoWithdata    <- GenoFreqs[!is.na(ActivityNumSum),]

Ethnicities     <- unique(GenoFreqs$Ethnicity)

# Get frequencies per activity level, locus and ethnicity
FreqData <- data.frame()
LocusCols <- c("Al512428", "Ac02980", "Ac021017")
Col <- "Ac02980"
for (Col in LocusCols){
  Loc <- GenoFreqs[,Col]
  
  SplitGeno <- sapply(Loc, function(x) strsplit(x, "/")[[1]])
  for (i in 1:nrow(SplitGeno)){
    SplitGeno[i, grep("-", SplitGeno[i,])] <- "0"
    SplitGeno[i,] <- gsub("*", "", SplitGeno[i,])
  }
  Activity  <- as.numeric(as.vector(t(SplitGeno)))
  Ethnicity <- rep(GenoFreqs$Ethnicity, 2)
  ConTab    <- table(Activity, Ethnicity)
  NewData   <- as.data.frame(ConTab / colSums(ConTab))
  NewData$Locus <- Col
  FreqData <- rbind(FreqData, NewData)
}

# Turn activities numeric
ActivityNum <- as.numeric(as.character(FreqData$Activity))

# Function to create samples
CreateSamples <- function(Ethnicity, Gender = "male", SampleSize){
  if (Gender == "male"){
    Loci <- rep(LocusCols,c(2, 1, 2))
  } else {
    Loci <- rep(LocusCols, 2)
  }
  blnEthnicity <- FreqData$Ethnicity == Ethnicity
  
  # Sample activity by locus
  ActByLocus <- sapply(Loci, function(x){
    blnSubset   <- (FreqData$Locus == x) & blnEthnicity 
    Probs       <- FreqData$Freq[blnSubset]
    Activities  <- ActivityNum[blnSubset]
    sample(Activities, SampleSize, replace = T,  prob = Probs)
  })
  
  # Return activity sums
  rowSums(ActByLocus)
}

# Create samples of activity sums
# SampleSizeBase  <- 1000
# Samples_males   <- c()
# Samples_females <- c()
# cat("Sampling L1 activity sums\n")
# for(x in Ethnicities) {
#   NrMale   <- sum(!GenoWithdata$isFemale & GenoWithdata$Ethnicity == x)
#   NrFemale <- sum(GenoWithdata$isFemale  & GenoWithdata$Ethnicity == x)
#   Samples_males <- c(Samples_males, 
#                      CreateSamples(x, "male", NrMale * SampleSizeBase))
#   Samples_females <- c(Samples_females, 
#                        CreateSamples(x, "female", NrFemale * SampleSizeBase))
# }

# Create qq plots for males and females
ActSum_males   <- ActivityNumSum[!is.na(ActivityNumSum) & !GenoFreqs$isFemale]
ActSum_females <- ActivityNumSum[GenoFreqs$isFemale]

# qqplot(Samples_males, ActSum_males, xlab = "Predicted", ylab = "Observed")
# lines(c(0, 1000), c(0, 1000))
# qqplot(Samples_females, ActSum_females, xlab = "Predicted", ylab = "Observed")
# lines(c(0, 1000), c(0, 1000))
# qqplot(c(Samples_females, Samples_males), c(ActSum_females, ActSum_males),
#        xlab = "Predicted", ylab = "Observed")
# lines(c(0, 1000), c(0, 1000))
# 
# max(Samples_females)
# hist(ActSum_males, breaks = seq(0, 400, 25))
# hist(ActSum_females, breaks = seq(0, 400, 25))
# hist(c(ActSum_females, ActSum_males), breaks = seq(0, 400, 25))
# hist(c(Samples_females, Samples_males), breaks = seq(0, 650, 25))

# Create 1000 samples of the size of the data
SampleSizeBase <- 1
cat("Sampling L1 activity sums with population structure\n")
StructuredSamples <- sapply(1:1000, function(x){
  Samples   <- c()
  for(x in Ethnicities) {
    NrMale   <- sum(!GenoWithdata$isFemale & GenoWithdata$Ethnicity == x)
    NrFemale <- sum(GenoWithdata$isFemale  & GenoWithdata$Ethnicity == x)
    SmallSamples_males   <- CreateSamples(x, "male", NrMale * SampleSizeBase)
    SmallSamples_females <- CreateSamples(x, "female", NrFemale * SampleSizeBase)
    Samples <- c(Samples, SmallSamples_males, SmallSamples_females)
  }
  Samples
})

# Function to calculate squared deviation between the cumulative distribution
# function of two samples
Sample1 <- c(ActSum_females, ActSum_males)
Sample2 <- c(Samples_females, Samples_males)
FsquaredDiff <- function(Sample1, Sample2, xVals = NULL,
                         NrIntervals = 1000){
  if (is.null(xVals)) {
    xMin <- min(c(Sample1, Sample2))
    xMax <- max(c(Sample1, Sample2))
    xVals <- seq(xMin, xMax, (xMax - xMin) / NrIntervals)
  }
  SquaredDiff <- sapply(xVals, function(x){
    F1 <- sum(Sample1 <= x) / length(Sample1)
    F2 <- sum(Sample2 <= x) / length(Sample2)
    (F1 - F2)^2
  })
  sum(SquaredDiff)
}

# Calculate squared differences between observations and sampled data
AllSamples <- c(Samples_females, Samples_males)
xVals <- 0:610
SqDiffActual <- FsquaredDiff(c(ActSum_females, ActSum_males),
                             AllSamples, xVals)
SqDiffSample <- sapply(1:ncol(StructuredSamples),function(x){
  FsquaredDiff(StructuredSamples[,x], AllSamples, xVals)
})
hist(SqDiffSample)
sum(SqDiffSample > SqDiffActual) / length(SqDiffSample)


# Calculate the variance of each sample and compare it to the actual variance
ActualVar     <- var(c(ActSum_males, ActSum_females))
ActualMean    <- mean(c(ActSum_males, ActSum_females))
SampledVars  <- apply(StructuredSamples, 2, var)
SampledMeans <- apply(StructuredSamples, 2, mean)
sum(SampledVars >= ActualVar) / length(SampledVars)
mean(SampledVars)
mean(SampledMeans)
ActualMean
dim(StructuredSamples)
