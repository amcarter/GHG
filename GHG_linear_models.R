# GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(MASS)
library(HH)
library(MuMIn)
library(lme4)
library(nlme)
library(lmerTest)
library(car)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")

# load slopes from whitebox
# slope <- read_csv('data/sites_whitebox_slopes.csv') %>%
#   slice(c(1:4,6,7)) %>%
#   mutate(site = c("NHC", 'PM', 'CBP', 'WB', 'WBP','UNHC')) %>%
#   left_join(dat, by = 'site') %>%
#   dplyr::select(site, slope_wbx = slope, slope_nhd) %>%
#   group_by(site) %>%
#   summarize_all(mean, na.rm = T) %>%
#   # whitebox slopes need to be rescaled because they were calculated on grouped cells (of 6)
#   # so the degrees are off in the original measurements:
#   mutate(slope_deg = (atan(6 * tan(slope_wbx*pi/180)) * (180/pi)),
#          slope_mm = tan(slope_deg * pi/180))
#
# write_csv(slope, 'data/sites_slope_comparison.csv')
slope <- read_csv('data/sites_slope_comparison.csv')

dat <- dat %>%
  left_join(slope[,c(1:2)])

ggplot(slope, aes(slope_nhd, slope_mm, col = site)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

png('figures/gas_evasion_coef_by_site.png', width = 6, height = 3, res = 300,
    units = 'in', family = 'cairo')
filter(dat, site !='MC751') %>%
  mutate(date = as.Date(group),
         site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC'))) %>%
ggplot(aes(site, K600)) +
  geom_boxplot(fill = 'grey') +
  ylab('K600 (day-1)')+
  ggtitle('Gas evasion coefficients by site') +
  theme_bw()
dev.off()
# PCA ####
# d2 <- mutate(d2, group = factor(group))
# d.inst <- d2 %>% select(site, group, ends_with(c('ugL', 'ugld', 'inst')),
#                         slope_nhd, habitat, GPP, ER, K600)
# d.avg <- d2 %>% select(site, group, habitat, ends_with(c('ugL', 'ugld')), depth,
#                        DO.obs, DO.sat, watertemp_C, discharge, slope_nhd,
#                        GPP, ER, K600)
#
# dat.pca.inst <- d.inst %>%
#   select(-site, -group, -habitat, -K600, -ends_with("ugld")) %>%
#   prcomp()
#
# autoplot(dat.pca.inst, data = d.inst, colour = "site", size = 2)
#
# dat.pca.avg <- d.avg %>%
#   select(-site, -group, -habitat, -depth, -ends_with("ugld")) %>%
#   prcomp()
#
# autoplot(dat.pca.avg, data = d.avg, colour = "discharge", size = 2)
#

# linear mixed effects models ####

# rescale covariates, normalize each to the mean
dat$CO2.flux_mgm2d <- dat$CO2.flux_ugld*dat$depth
dat$CH4.flux_mgm2d <- dat$CH4.flux_ugld*dat$depth
dat$O2.flux_mgm2d <- dat$O2.flux_ugld*dat$depth
dat$N2O.flux_mgm2d <- dat$N2O.flux_ugld*dat$depth
scaled <- dat %>%
    filter(site != 'MC751',
         !is.na(CH4.ugL),
         !is.na(GPP)) %>%
    mutate(logQ = log(discharge),
         NER = ER - GPP,
         DO.persat = DO.obs/DO.sat,
         # logWRT = log(1000*depth*width_march_m/discharge),
         no3n_mgl = ifelse(no3n_mgl == 0, 0.0015, no3n_mgl), # replace zero with mdl
         site = factor(site, levels=c('UNHC','WBP','WB','CBP','PM','NHC')),
         log_no3n = log(no3n_mgl)
         ) %>%
    dplyr::select(site, logQ, watertemp_C, ER, GPP, DO.obs, slope_deg, depth,
            log_no3n, doc_mgl, ends_with(c("ugld", "ugL", 'mgm2d'))) %>%
    mutate(across(-c(site), ~ scale(.)[,1, drop = T]),
         across(c(site), ~ factor(.)))


preds <- scaled %>%
  dplyr::select(logQ, watertemp_C, GPP, ER, slope_deg, DO.obs,
         log_no3n, doc_mgl, depth)

pred_cov <- data.frame(cov(preds))
write_csv(pred_cov, 'data/linear_models/predictor_covariance_matrix.csv')

cor.test(preds$doc_mgl, preds$depth, method = 'pearson')

PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)

  return(PRESS)
}

pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)

  return(pred.r.squared)
}


#testing if random effects are important ####
  #is the correlation within sites higher than the correlation between sites?

names(scaled)
#random intercepts model
lme0site <- lme(CH4.ugL~1, random=~1|site, data=scaled)
# lme0site <- lme4::lmer(CH4.ugL~1 + (1|site), data=scaled)
summary(lme0site)
VarCorr(lme0site)
tau.sq<-as.numeric(VarCorr(lme0site)[1,1])
sigma.sq<-as.numeric(VarCorr(lme0site)[2,1])
tau.sq/(tau.sq+sigma.sq)
#correlation between observations at the same site is only 0.01
#probably doesn't justify a lme model - at least that makes it easier!
# test if a random effects model is justified
#one thought was that GHG patterns by site were more stable than between sites
#so should probably test the site as a fixed effect

summary(aov(CO2.ugL ~ site, data = scaled))         # not significant
summary(aov(CO2.flux_mgm2d ~ site, data = scaled))   # not significant
summary(aov(CH4.ugL ~ site, data = scaled))         # ** p = 0.00324
summary(aov(CH4.flux_mgm2d ~ site, data = scaled))   # p = 0.0687
summary(aov(N2O.ugL ~ site, data = scaled))         # not significant
summary(aov(N2O.flux_mgm2d ~ site, data = scaled))   # * p = 0.0347
summary(aov(DO.obs ~ site, data = scaled))          # not significant
summary(aov(O2.flux_mgm2d ~ site, data = scaled))    # not significant

summary(lm(CH4.ugL ~ site-1, data = scaled))
summary(lm(CH4.flux_ugld ~ site-1, data = scaled))
summary(lm(N2O.flux_ugld ~ site-1, data = scaled))

# Use Leaps package to find the best models for each gas concentration and flux
find_best_model <- function(preds, y, gas, flux = FALSE){
    if(gas != 'N2O') preds <- dplyr::select(preds, -log_no3n)
    if(gas == 'O2') preds <- dplyr::select(preds, -DO.obs)
    if(flux) preds <- dplyr::select(preds, -depth)
    a <- leaps::regsubsets(x = as.matrix(preds), y = y, nbest = 8, nvmax = 5)
    asum <- summary(a)

    mods <- bind_cols(asum$which, rsq = asum$rsq, adjr2 = asum$adjr2,
                      cp = asum$cp, bic = asum$bic) %>%
        dplyr::select(-'(Intercept)')
    w <- which(mods$logQ & mods$doc_mgl)
    if(isFALSE(flux)) w <- c(w, which(mods$depth & mods$doc_mgl))
    if(gas != "O2") w <- c(w, which(mods$DO.obs & mods$ER))

    mods <- mods[-w,] %>%
        arrange(bic, cp, -adjr2) %>%
        mutate(aicc = NA_real_,
               rsq_pred = NA_real_)
    mods$model = seq(1:nrow(mods))

    preds$y <- y
    coefs <- data.frame()
    for(i in 1:nrow(mods)){
        vars <- mods[i,] %>%
            dplyr::select(-rsq, -adjr2,-cp,-bic, -aicc, -rsq_pred) %>%
            dplyr::select_if(isTRUE) %>%
            colnames()
        m <- lm(paste0('y ~ ',
                       paste(vars, collapse=' + ')),
                data = preds)
        mods$aicc[i] <- AICc(m)
        mods$rsq_pred[i] <- pred_r_squared(m)

        cc <- data.frame(summary(m)$coefficients) %>%
            slice(-1)
        cc <- mutate(cc, pred = rownames(cc)) %>%
            dplyr::select(pred, mean = Estimate,
                          se = 'Std..Error', p = 'Pr...t..') %>%
            mutate(model = i,
                   gas = gas,
                   flux = flux)
        row.names(cc) <- NULL

        coefs <- bind_rows(coefs, cc)
    }
    w <-  which(mods$aicc <= arrange(mods, aicc)$aicc[5])
    mods <- mods[w,] %>%
        arrange(aicc, -rsq_pred,-adjr2)
    coefs <- filter(coefs, model %in% unique(mods$model))

    mods <- mutate(mods,
                   delta_aicc = aicc - min(aicc),
                   delta_rsqp = rsq_pred - max(rsq_pred),
                   rel_likelihood = exp(-0.5 * delta_aicc),
                   aic_weight = rel_likelihood/sum(rel_likelihood),
                   gas = gas,
                   flux = flux)
    return(list(mods, coefs))
}

CO2.mods <- find_best_model(preds, scaled$CO2.ugL, 'CO2')
CO2flux.mods <- find_best_model(preds, scaled$CO2.flux_mgm2d, 'CO2', TRUE)
CH4.mods <- find_best_model(preds, scaled$CH4.ugL, 'CH4')
CH4flux.mods <- find_best_model(preds, scaled$CH4.flux_mgm2d, 'CH4', TRUE)
N2O.mods <- find_best_model(preds, scaled$N2O.ugL, 'N2O')
N2Oflux.mods <- find_best_model(preds, scaled$N2O.flux_mgm2d, 'N2O', TRUE)
O2.mods <- find_best_model(preds, scaled$DO.obs, 'O2')
O2flux.mods <- find_best_model(preds, scaled$O2.flux_mgm2d, 'O2', TRUE)

best_mods <- bind_rows(N2O.mods[[1]], N2Oflux.mods[[1]],
                       CO2.mods[[1]], CO2flux.mods[[1]],
                       CH4.mods[[1]], CH4flux.mods[[1]],
                       O2.mods[[1]], O2flux.mods[[1]])
mod_coeffs <- bind_rows(N2O.mods[[2]], N2Oflux.mods[[2]],
                        CO2.mods[[2]], CO2flux.mods[[2]],
                        CH4.mods[[2]], CH4flux.mods[[2]],
                        O2.mods[[2]], O2flux.mods[[2]])

best_mod_coeffs <- best_mods[c(1,6,11,16,21,26),] %>%
    dplyr::select(model, gas, flux) %>%
    left_join(mod_coeffs) %>%
    dplyr::select(-model)

write_csv(best_mods, 'data/linear_models/best_linear_model_summaries.csv')
write_csv(best_mod_coeffs, 'data/linear_models/best_linear_model_coefficients.csv')
