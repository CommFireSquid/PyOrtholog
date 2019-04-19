setwd('D:/Data')
CTR <- read.csv('CTR.csv')
GLY <- read.csv('GLY.csv')
ELM <- read.csv('ELM.csv')
SLM <- read.csv('SLM.csv')
FAB <- read.csv('FAB.csv')
TRY <- read.csv('TRY.csv')

# Test population Variance and Population Mean Individually VS Control

# 1st merge the dataframes
CTR_GLY = merge(CTR,GLY,by="Species")
CTR_ELM = merge(CTR,ELM,by="Species")
CTR_SLM = merge(CTR,SLM,by="Species")
CTR_FAB = merge(CTR,FAB,by="Species")
CTR_TRY = merge(CTR,TRY,by="Species")

# Hypothesis: There is a significant variation between each set of paired measures
t.test(CTR_GLY$Average.x, CTR_GLY$Average.y, alternative = "two.sided", paired=TRUE)
t.test(CTR_ELM$Average.x, CTR_ELM$Average.y, alternative = "two.sided", paired=TRUE)
t.test(CTR_SLM$Average.x, CTR_SLM$Average.y, alternative = "two.sided", paired=TRUE)
t.test(CTR_FAB$Average.x, CTR_FAB$Average.y, alternative = "two.sided", paired=TRUE)
t.test(CTR_TRY$Average.x, CTR_TRY$Average.y, alternative = "two.sided", paired=TRUE)

# Hypothesis: In all cases, the pathway relationship indicates a lower genetic distance
t.test(CTR_GLY$Average.y, CTR_GLY$Average.x, alternative = "less", paired=TRUE)
t.test(CTR_ELM$Average.y, CTR_ELM$Average.x, alternative = "less", paired=TRUE)
t.test(CTR_SLM$Average.y, CTR_SLM$Average.x, alternative = "less", paired=TRUE)
t.test(CTR_FAB$Average.y, CTR_FAB$Average.x, alternative = "less", paired=TRUE)
t.test(CTR_TRY$Average.y, CTR_TRY$Average.x, alternative = "less", paired=TRUE)