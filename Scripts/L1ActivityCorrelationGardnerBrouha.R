# Read in table by Gardener et al. 2017
GardenerTable <- read.csv("D:/OneDrive - American University of Beirut/L1polymORF/Data/Gardner et al. 2017 Supp Table 9.csv",
                          as.is = T)

plot(GardenerTable$Total.Offspring, GardenerTable$Activity.of.L1RP..Brouha.2002..2003..or.L1.3..Beck.2010.)
cor.test(GardenerTable$Total.Offspring, GardenerTable$Activity.of.L1RP..Brouha.2002..2003..or.L1.3..Beck.2010.)
cor.test(GardenerTable$Total.Offspring, 
         GardenerTable$Activity.of.L1RP..Brouha.2002..2003..or.L1.3..Beck.2010.,
         method = "spearman")
