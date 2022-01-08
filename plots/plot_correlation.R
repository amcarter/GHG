# Make correlation plot for GHG and drivers

library(tidyverse)
library(corrplot)
setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
vars <- dat %>%
  select(ends_with('.ugl'), O2.mgL = DO.obs, 
         watertemp = watertemp_C, slope = slope_nhd, GPP, ER, 
         NO3.mgl = no3n_mgl, DOC.mgl = doc_mgl, discharge)

res <- cor(vars, use='complete.obs')
res1 <- cor.mtest(vars, conf.level=0.95)

write.csv(res, 'data/GHG_correlationtable.csv')
write.csv(res1, 'data/corsignificance.csv')

pal1<-colorRampPalette(c('coral2', 'white', 'steelblue'))

tiff('figures/correlationplot.tif', width=6, height=6, res=390, units='in')
corrplot(res, type = "lower", order = "original", p.mat = res1$p,
         tl.col = "black", tl.srt = 90, method='color', diag=FALSE,
         #  addCoef.col = "black",
         col=pal1(20),
         sig.level = c(0.01, 0.05, 0.1), insig = 'label_sig',
         pch.cex = 1, pch.col='grey20')
mtext('Correlation significance', adj = .79)
mtext(expression('* ' <= ' 0.1 '), adj = .8055, line = -1)
mtext(expression(' ** ' <= ' 0.05'), adj = .808, line = -2)
mtext(expression('*** ' <= ' 0.01'), adj = .808, line = -3)
dev.off()

