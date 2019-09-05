# The following script explores the correlation between L1 frequency and 
# retrotransposition activity
L1Catalog <- read.csv("D:/L1polymORF/Data/L1CatalogExtended.csv")
cor.test(L1Catalog$ActivityNum, L1Catalog$Allele_frequency_Num, 
         use = "pairwise.complete.obs", method = "spearman")
plot(L1Catalog$ActivityNum, L1Catalog$Allele_frequency_Num)
blnFreq <- L1Catalog$Allele_frequency_Num < 1
cor.test(L1Catalog$ActivityNum[blnFreq], L1Catalog$Allele_frequency_Num[blnFreq], 
         use = "pairwise.complete.obs", method = "spearman")
cor.test(L1Catalog$ActivityNum[blnFreq], L1Catalog$Allele_frequency_Num[blnFreq], 
         use = "pairwise.complete.obs")
