# The following script goes through sampled model fits and calculates the mean likelihood
# average over all sampled data
load("D:/OneDrive - American University of Beirut/L1polymORF/Data/L1SelectionResults_MELT_GroupwithSim.RData")

ModelLogLiks <- sapply(ModelFitSamples_pracma, function(x){
  c(LLa  = x$ML_a$value,
    LLab  = x$ML_ab$value,
    LLac  = x$ML_ac$value,
    LLabc = x$ML_abc$value)
})
rowMeans(ModelLogLiks)

ModelAics <- sapply(ModelFitSamples_pracma, function(x){
  as.numeric(as.character(x$AICTab$AIC))
})
meanAIC <- rowMeans(ModelAics)
meanAIC - min(meanAIC)
HistTrue <- hist(L1DetectAgg_withL1$L1widthTrue, plot = F)
HistEst  <- hist(L1DetectAgg_withL1$L1widthEst, plot = F)
plot(HistTrue$density, HistEst$density)
text(HistTrue$density, HistEst$density, HistTrue$mids)
lines(c(0, 1), c(0, 1))
cor(HistTrue$density, HistEst$density)
HistEst
table()