# Make correlation plot for GHG and drivers

library(tidyverse)
library(corrplot)
setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
vars <- dat %>%
  select(ends_with('.ugl'), O2.mgL = DO.obs,
         watertemp = watertemp_C, slope = slope_mm, depth, GPP, ER,
         NO3.mgl = no3n_mgl, DOC.mgl = doc_mgl, discharge)
fluxvars <- dat %>%
  select(ends_with('_ugld'),
         watertemp = watertemp_C, slope = slope_mm, depth, GPP, ER,
         NO3.mgl = no3n_mgl, DOC.mgl = doc_mgl, discharge)

res <- cor(vars, use='complete.obs')
res1 <- cor.mtest(vars, conf.level=0.95)
fluxres <- cor(fluxvars, use='complete.obs')
fluxres1 <- cor.mtest(fluxvars, conf.level=0.95)


dd <- vars %>%
    mutate(NO3.mgl = ifelse(NO3.mgl == 0, 0.0015, NO3.mgl),
           across(any_of(c('discharge', 'DOC.mgl', 'NO3.mgl')), ~log(.))) %>%
    select(slope, depth, 'log(Q)' = discharge, 'water temp' = watertemp,
           GPP, ER, ':O[2]' = O2.mgL, 'log(DOC)' = DOC.mgl,
           ':log(NO[3])' = NO3.mgl)
d <- cor(dd, use = 'complete.obs')
d1 <- cor.mtest(dd, conf.level = 0.95)

write.csv(res, 'data/GHG_correlationtable.csv')
write.csv(res1, 'data/corsignificance.csv')

pal1<-colorRampPalette(c('coral2', 'white', 'steelblue'))

tiff('figures/final/correlationplot.tif', width = 4, height = 4,
    res = 300, units = 'in')
    par(oma = c(.5,0,0,0))
    corrplot(d, type = "lower", order = "original", p.mat = d1$p,
             tl.col = "black", tl.srt = 90, tl.cex = .7,
             method='color', diag=FALSE, cl.cex = .6,
             cl.length = 11, col=pal1(20),
             sig.level = c(0.01, 0.05, 0.1), insig = 'label_sig',
             pch.cex = 1, pch.col='grey20')
    mtext('Correlation Coefficient', 1, 4.5, cex = .7, adj = .57)

dev.off()

tiff('figures/final/correlationplot_2panel.tif', width=12, height=6,
     res=390, units='in')
    par(mfrow = c(1,2))
    corrplot(res, type = "lower", order = "original", p.mat = res1$p,
             tl.col = "black", tl.srt = 90, method='color', diag=FALSE,
             #  addCoef.col = "black",
             col=pal1(20),
             sig.level = c(0.01, 0.05, 0.1), insig = 'label_sig',
             pch.cex = 1, pch.col='grey20')
    mtext('Gas Concentrations', 3, adj = 0.8)
    corrplot(fluxres, type = "lower", order = "original", p.mat = fluxres1$p,
             tl.col = "black", tl.srt = 90, method='color', diag=FALSE,
             #  addCoef.col = "black",
             col=pal1(20),
             sig.level = c(0.01, 0.05, 0.1), insig = 'label_sig',
             pch.cex = 1, pch.col='grey20')
    mtext('Gas Fluxes', 3, adj = 0.8)
    # mtext('Correlation significance', adj = .79, line = -2, cex = .8)
    # mtext(expression('* ' <= ' 0.1 '), adj = .8055, line = -3, cex = .8)
    # mtext(expression(' ** ' <= ' 0.05'), adj = .808, line = -4, cex = .8)
    # mtext(expression('*** ' <= ' 0.01'), adj = .808, line = -5,  cex = .8)
dev.off()

