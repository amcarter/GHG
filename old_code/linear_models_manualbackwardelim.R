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
# setwd('C://Users/adelv/Dropbox/Duke/NHC/Alice')


dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")

# load slopes from whitebox
# slope <- read_csv('data/sites_whitebox_slopes.csv') %>%
#   slice(c(1:4,6,7)) %>%
#   mutate(site = c("NHC", 'PM', 'CBP', 'WB', 'WBP','UNHC')) %>%
#   left_join(dat, by = 'site') %>%
#   dplyr::select(site, slope_wbx = slope, slope_nhd) %>%
#   group_by(site) %>%
#   summarize_all(mean, na.rm = T) %>%
#   mutate(slope_2deg = (atan(6 * tan(slope_wbx*pi/180)) * (180/pi)),
#          slope = tan(slope_2deg * pi/180))
#
# write_csv(slope, 'data/sites_slope_comparison.csv')
slope <- read_csv('data/sites_slope_comparison.csv')

dat <- dat %>%
  left_join(slope[,c(1:2)]) %>%
  rename(slope_deg = slope_wbx)

ggplot(slope, aes(slope_nhd, slope, col = site)) + geom_point() +
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

get_mod_results <- function(m, gas, flux = F){
  # m <- stepAIC(mm, direction = 'both', trace = F)
  # print("Variance inflation factors - should be less than ~5")
  # print(vif(m))
  r <- data.frame(r.squaredGLMM(m)) %>%
    mutate(gas = gas, flux = flux,
           formula = as.character(summary(m)$call)[2])
  r$pred.rsq = pred_r_squared(m)
  cf <- summary(m)$coefficients %>%
    as_tibble() %>%
    mutate(predictor = rownames(summary(m)$coefficients))
  colnames(cf) <- c('estimate', 'std_error', 't_val',
                    'pr_greaterthan_t', 'predictor')
  ci <- data.frame(confint(m)) %>%
    mutate(predictor = rownames(confint(m))) %>%
    left_join(cf, by = 'predictor') %>%
    mutate(effect = "fixed",
           gas = gas,
           flux = flux)

  # steps <- m$anova %>%
  #   as_tibble() %>%
  #   rename(Resid_Df = 'Resid. Df', resid_dev = 'Resid. Dev') %>%
  #   mutate(effect = "fixed",
  #          gas = gas,
  #          flux = flux)


    return(list(r2 = r,
                # steps = steps,
                mod = ci))
}

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

# manual model stepping ####
#I like to step through models manually because sometimes a really small drop
#in the AIC involves removing a variable close to significant
mod_steps <- data.frame()

# CO2.conc ####
mm <- lm(CO2.ugL ~ site + slope_deg + logQ + watertemp_C +
           GPP + ER + DO.obs + doc_mgl, data = scaled)
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
summary(mm) # 0.835
pred_r_squared(mm)# .723
AICc(mm) #75.08
mm <- lm(CO2.ugL ~ slope_deg + logQ + watertemp_C +
           GPP + ER + DO.obs + doc_mgl, data = scaled)
summary(mm) # 0.816
pred_r_squared(mm)# .731
AICc(mm) #71.92
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.ugL ~ slope_deg + logQ + watertemp_C +
           GPP + ER + DO.obs , data = scaled)
summary(mm) # 0.82
pred_r_squared(mm)# .742
AICc(mm) #69.1
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.ugL ~ slope_deg + logQ + watertemp_C +
           GPP + DO.obs , data = scaled)
summary(mm) # 0.824
pred_r_squared(mm)# .775
AICc(mm) #66.3
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.ugL ~ slope_deg + watertemp_C + GPP + DO.obs , data = scaled)
summary(mm) # 0.819
pred_r_squared(mm)# .774
AICc(mm) #66.1
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.ugL ~ slope_deg + watertemp_C + DO.obs , data = scaled)
summary(mm) # 0.813
pred_r_squared(mm)# .782
AICc(mm) #66.3
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.ugL ~ watertemp_C + DO.obs , data = scaled)
summary(mm) # 0.811
pred_r_squared(mm)# .785
AICc(mm) #65.4
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
CO2.mod <- lm(CO2.ugL ~  watertemp_C + DO.obs, data = scaled)

mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "CO2", unit = "conc") %>%
  bind_rows(mod_steps)

mm <- lm(CO2.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl, data = scaled)

summary(mm) #.795
AICc(mm) #85.9
pred_r_squared(mm)#.725
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs, data = scaled)
summary(mm) #.789
AICc(mm) #85.2
pred_r_squared(mm)#.725
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER, data = scaled)
summary(mm) #.793
AICc(mm) #82.1
pred_r_squared(mm)#.728
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.flux_mgm2d ~ site + logQ + watertemp_C +
             GPP + ER, data = scaled)
summary(mm) #.793
AICc(mm) #82.1
pred_r_squared(mm)#.728
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.flux_mgm2d ~ site + logQ + GPP + ER, data = scaled)
summary(mm) #.787
AICc(mm) #81.3
pred_r_squared(mm)#.726
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CO2.flux_mgm2d ~ logQ + GPP + ER, data = scaled)
summary(mm) #.745
AICc(mm) #82.1
pred_r_squared(mm)#.712
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
CO2flux.mod <- lm(CO2.flux_mgm2d ~ logQ + GPP + ER , data = scaled)

mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "CO2", unit = "flux") %>%
  bind_rows(mod_steps)

mm <- lm(CH4.ugL ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl, data = scaled)
summary(mm) #.775
AICc(mm) #90.8
pred_r_squared(mm)#.607
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.ugL ~ slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl, data = scaled)
summary(mm) #.76
AICc(mm) #85.5
pred_r_squared(mm)#.607
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.ugL ~ slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs , data = scaled)
summary(mm) #.76
AICc(mm) #83.7
pred_r_squared(mm)#.613
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.ugL ~ slope_deg + logQ + watertemp_C +
             GPP + DO.obs , data = scaled)
summary(mm) #.749
AICc(mm) #84.4
pred_r_squared(mm)#.643
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.ugL ~ slope_deg + logQ + GPP + DO.obs , data = scaled)
summary(mm) #.748
AICc(mm) #83.1
pred_r_squared(mm)#.657
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
CH4.mod <- lm(CH4.ugL ~ slope_deg + logQ + GPP + DO.obs , data = scaled)

mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "CH4", unit = "conc") %>%
  bind_rows(mod_steps)


mm <- lm(CH4.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl, data = scaled)
summary(mm) #.667
AICc(mm) #110.7
pred_r_squared(mm)#.554
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.flux_mgm2d ~ site + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl, data = scaled)
summary(mm) #.667
AICc(mm) #110.7
pred_r_squared(mm)#.554
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.flux_mgm2d ~ site + logQ + watertemp_C +
             GPP + ER + DO.obs , data = scaled)
summary(mm) #.675
AICc(mm) #107
pred_r_squared(mm)#.582
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.flux_mgm2d ~ site + logQ + watertemp_C + GPP + ER, data = scaled)
summary(mm) #.681
AICc(mm) #104
pred_r_squared(mm)#.601
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.flux_mgm2d ~ site + logQ + GPP + ER, data = scaled)
summary(mm) #.675
AICc(mm) #103
pred_r_squared(mm)#.596
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(CH4.flux_mgm2d ~ site  + GPP + ER, data = scaled)
summary(mm) #.657
AICc(mm) #103.8
pred_r_squared(mm)#.574
mod <-data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))

CH4flux.mod <- lm(CH4.flux_mgm2d ~ site + GPP + ER, data = scaled)
mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "CH4", unit = "flux") %>%
  bind_rows(mod_steps)


mm <- lm(N2O.ugL ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.569
AICc(mm) #126
pred_r_squared(mm)#.392
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.ugL ~ slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.587
AICc(mm) #115.2
pred_r_squared(mm)#.478
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.ugL ~ slope_deg + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.596
AICc(mm) #112.2
pred_r_squared(mm)#.507
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.ugL ~ slope_deg + watertemp_C +
             GPP + ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.574
AICc(mm) #113
pred_r_squared(mm)#.487
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.ugL ~ slope_deg + watertemp_C +
             ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.582
AICc(mm) #110.4
pred_r_squared(mm)#.515
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.ugL ~ slope_deg + ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.587
AICc(mm) #108.2
pred_r_squared(mm)#.528
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
N2O.mod <- lm(N2O.ugL ~ slope_deg + ER + doc_mgl + log_no3n, data = scaled)
mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "N2O", unit = "conc") %>%
  bind_rows(mod_steps)


mm <- lm(N2O.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.55
AICc(mm) #128
pred_r_squared(mm)#.392
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.flux_mgm2d ~ slope_deg + logQ + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.54
AICc(mm) #120
pred_r_squared(mm)#.45
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.flux_mgm2d ~ slope_deg  + watertemp_C +
             GPP + ER + DO.obs + doc_mgl + log_no3n, data = scaled)
summary(mm) #.548
AICc(mm) #117.9
pred_r_squared(mm)#.477
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.flux_mgm2d ~ slope_deg  + watertemp_C +
             GPP + ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.48
AICc(mm) #123
pred_r_squared(mm)#.407
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.flux_mgm2d ~ slope_deg  + watertemp_C +
             ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.493
AICc(mm) #120.3
pred_r_squared(mm)#.433
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(N2O.flux_mgm2d ~ slope_deg + ER + doc_mgl + log_no3n, data = scaled)
summary(mm) #.483
AICc(mm) #116.7
pred_r_squared(mm)#.423
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
N2Oflux.mod <- lm(N2O.flux_mgm2d ~ slope_deg + ER + doc_mgl + log_no3n, data = scaled)
mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "N2O", unit = "flux") %>%
  bind_rows(mod_steps)


mm <- lm(DO.obs ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + doc_mgl, data = scaled)
summary(mm) #.56
AICc(mm) #122.38
pred_r_squared(mm)#.384
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(DO.obs ~  slope_deg + logQ + watertemp_C +
             GPP + ER + doc_mgl, data = scaled)
summary(mm) #.579
AICc(mm) #112.5
pred_r_squared(mm)#.452
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(DO.obs ~  slope_deg + logQ + watertemp_C +
             GPP + ER, data = scaled)
summary(mm) #.548
AICc(mm) #114.4
pred_r_squared(mm)#.45
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(DO.obs ~  slope_deg + watertemp_C +
             GPP + ER, data = scaled)
summary(mm) #.525
AICc(mm) #115.3
pred_r_squared(mm) #.44
mod <- data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
O2.mod <- lm(DO.obs ~ slope_deg + watertemp_C + GPP + ER, data = scaled)
mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "O2", unit = "conc") %>%
  bind_rows(mod_steps)


mm <- lm(O2.flux_mgm2d ~ site + slope_deg + logQ + watertemp_C +
             GPP + ER + doc_mgl, data = scaled)
summary(mm) #.873
AICc(mm) #59
pred_r_squared(mm)#.807
mod <-  data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(O2.flux_mgm2d ~ slope_deg + logQ + watertemp_C +
             GPP + ER + doc_mgl, data = scaled)
summary(mm) #.883
AICc(mm) #47.3
pred_r_squared(mm)#.837
mod <-  data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(O2.flux_mgm2d ~ slope_deg + logQ + watertemp_C +
             GPP + ER , data = scaled)
summary(mm) #.885
AICc(mm) #44.6
pred_r_squared(mm)#.84
mod <-  data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(O2.flux_mgm2d ~ logQ + watertemp_C + GPP + ER , data = scaled)
summary(mm) #.888
AICc(mm) #41.9
pred_r_squared(mm)#.845
mod <-  data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
mm <- lm(O2.flux_mgm2d ~  watertemp_C + GPP + ER , data = scaled)
summary(mm) #.881
AICc(mm) #43.1
pred_r_squared(mm)#.838
mod <-  data.frame(model = as.character(mm$call)[2],
                            adj_r2 = summary(mm)$adj.r.squared,
                            pred_r2 = pred_r_squared(mm),
                            AICc = AICc(mm))
O2flux.mod <- lm(O2.flux_mgm2d ~  watertemp_C + GPP + ER , data = scaled)
mod_steps <- mod %>%
  mutate(pr_adj_ratio = pred_r2/adj_r2,
         gas = "O2", unit = "flux") %>%
  bind_rows(mod_steps)

# compile backward selection ####
rsqs <- data.frame()
mod_fits <- data.frame()

handmods <- list(O2.mod, CO2.mod, CH4.mod, N2O.mod,
                 O2flux.mod, CO2flux.mod, CH4flux.mod, N2Oflux.mod)

gasses <- c('N2O', 'O2', 'CO2', 'CH4')
for(i in 1:8){

  if(i <= 4){ flux = FALSE} else { flux = TRUE}
  gas = gasses[(i %% 4) + 1]
  mm <- handmods[[i]]
  out <- get_mod_results(mm,gas = gas, flux = flux)
  out$r2$AIC <- AICc(mm)
  out$r2$adj_r2 <- summary(mm)$adj.r.squared
  print(vif(mm))
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)
}

saveRDS(list(mod_fits = mod_fits,
             rsqs = rsqs,
             mod_steps = mod_steps),
        'data/linear_models/model_fits_steps_rsqs_lm_manual_AIC_backward_elimination_03_05.rds')

write_csv(mod_steps, 'data/linear_models/backward_elimination_steps_0305.csv')
mod_fits %>%
  as_tibble() %>%
  mutate(predictor = ifelse(predictor == '(Intercept)', 'Intercept', predictor),
         flux = ifelse(flux, 'flux','conc'))%>%
  rename(p_val = pr_greaterthan_t, ci_2.5 = X2.5.., ci_97.5 = X97.5..)%>%
  dplyr::select(-effect) %>%
  pivot_wider(names_from = c(gas, flux),
              values_from = c(estimate, std_error,
                              ci_2.5, ci_97.5, t_val, p_val))%>%
write_csv('data/linear_models/backward_elimination_model_fits_0305.csv')
write_csv(rsqs, 'data/linear_models/backward_elimination_rsqs_0305.csv')
#CO2.conc ####
#remove discharge
mmQ <- lm(CO2.ugL ~ watertemp_C + depth +
           GPP + ER + DO.obs + doc_mgl + site, data = scaled)
Anova(mmQ)
summary(mmQ) #0.878
AIC(mmQ) #51.98, an improvement
vif(mmQ) #I haven't used this before and got an error - not  sure how to use this properly
#try removing site next

mmQ.site <- lm(CO2.ugL ~ watertemp_C + depth + GPP + ER + DO.obs + doc_mgl,
                data = scaled)
summary(mmQ.site) #.876
AIC(mmQ.site) #48.98 Better
vif(mmQ.site)

# try removing DO because it is correlated with ER
mmQ.site.DO <- lm(CO2.ugL ~ watertemp_C + depth + GPP + ER + doc_mgl,
                data = scaled)

AIC(mmQ.site.DO) #112 way worse
summary(mmQ.site.DO) #0.59 way worse

# everything else is significant, but try removing ER, the least sig
mmQ.site.ER <- lm(CO2.ugL ~ watertemp_C + depth + GPP + DO.obs + doc_mgl,
                  data = scaled)
AIC(mmQ.site.ER) # 50.39 worse
summary(mmQ.site.ER) #0.871
# also worse, but just barely. I think leave ER out bc correlation with DO

# final model:
CO2.mod <- mmQ.site.ER
r.squaredGLMM(CO2.mod)
pred_r_squared(CO2.mod)
#so might advocate for removing logQ for that slight improvement to AIC
#because of the strong colinearity issue (but discussing that in the text)
#this might be a good time to run the vif - I don't see other candidates for
#predictor removal based on cov unless the vif is high

#if date wont be added in (which I don't *think* it needs to be since temp
#especially is here) then prob an important point is dicussing driver change
#over time - which I think you do already, just stressing it

#last point is the overfitting issue.  I haven't dealt with that in a while!
#doing some quick reading, calculating a predicted R squared might be helpful
#https://stackoverflow.com/questions/64224124/how-to-calculate-predicted-r-sq-in-r


# CO2 flux ####
mm <- lm(CO2.flux_ugld ~ logQ + watertemp_C +  depth + GPP + ER +
           DO.obs + doc_mgl + site, data = scaled)
AIC(mm) # 80.034
Anova(mm)
summary(mm)# rsq = .7979
# doc looks definitely not useful:
mmdoc <- lm(CO2.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER +
              DO.obs + site, data = scaled)
anova(mmdoc)
summary(mmdoc) #rsq = .8026, definitely better
AIC(mmdoc)# 78.05
# next worst predictor is depth
mmdoc.d <- lm(CO2.flux_ugld ~ logQ + watertemp_C  + GPP + ER + DO.obs +
              site, data = scaled)
anova(mmdoc.d)
summary(mmdoc.d) #rsq = .8072
AIC(mmdoc.d)# 76.06, improvement
#next DO
mmdoc.d.DO <- lm(CO2.flux_ugld ~ logQ + watertemp_C +  GPP + ER +
              site, data = scaled)
Anova(mmdoc.d.DO)
summary(mmdoc.d.DO) #rsq = .809
AIC(mmdoc.d.DO)# 74.7, Improvement
#temperature
mmdoc.d.DO.temp <- lm(CO2.flux_ugld ~ logQ + GPP + ER +
                       site, data = scaled)
Anova(mmdoc.d.DO.temp)
summary(mmdoc.d.DO.temp ) #rsq = .812
AIC(mmdoc.d.DO.temp )# 72.966, improvement
pred_r_squared(mmdoc.d.DO.temp) # 0.764
vif(mmdoc.d.DO.temp)

# try discharge
mmdoc.d.DO.temp.Q <- lm(CO2.flux_ugld ~ GPP + ER +
                       site, data = scaled)
summary(mmdoc.d.DO.temp.Q) #rsq = .803
AIC(mmdoc.d.DO.temp.Q)# 75.0, worse
pred_r_squared(mmdoc.d.DO.temp.Q) # 0.763

# everything left is significant and vif is good, leave Q in
CO2flux.mod <- mmdoc.d.DO.temp
r.squaredGLMM(CO2flux.mod)
pred_r_squared(CO2flux.mod)
mmdoc.d.DO.temp <- lm(CO2.flux_ugld ~ logQ + GPP + ER +
                       site-1, data = scaled)
summary(mmdoc.d.DO.temp ) #rsq = .812

# CH4. conc ####
mm <- lm(CH4.ugL ~ logQ + watertemp_C + depth + GPP + ER + DO.obs +
           doc_mgl + site, data = scaled)
# mm <- lme(CH4.ugL ~ logQ + watertemp_C + GPP + ER + DO.obs +
#             doc_mgl + CO2.ugL + slope_nhd, random = ~1|site, data = scaled)
# library(MuMIn)
# r.squaredGLMM(mm)
Anova(mm)
summary(mm)# .799
AIC(mm) #79.737
stepAIC(mm)

# doc not sig in either test
mmdoc <- lm(CH4.ugL ~ logQ + watertemp_C + depth + GPP + ER +
              DO.obs + site, data = scaled)
Anova(mmdoc)
summary(mmdoc)# .8037
AIC(mmdoc) # 77.74, improvement
vif(mmdoc)

# water temp
mmdoc.temp <- lm(CH4.ugL ~ logQ + depth + GPP + ER +
                 DO.obs + site, data = scaled)
anova(mmdoc.temp)
summary(mmdoc.temp)# .799
AIC(mmdoc.temp) # 78.35, worse
# remove depth
mmdoc.d <- lm(CH4.ugL ~ logQ + watertemp_C + GPP + ER + DO.obs +
                   site, data = scaled)
Anova(mmdoc.d)
summary(mmdoc.d)# .797
AIC(mmdoc.d) # 77.72, worse, but not by much. Does it improve the VIF issue?
vif(mmdoc.d) # yep
pred_r_squared(mmdoc.d) #.675
#try site
mmdoc.site <- lm(CH4.ugL ~ logQ + watertemp_C +depth + ER + GPP + DO.obs,
                           data = scaled)
summary(mmdoc.site)# .707
AIC(mmdoc.site) # 95.35 much worse

# try temp again, or combinations of temp + GPP, ER
mmdoc.d.temp <- lm(CH4.ugL ~ logQ + GPP + ER + DO.obs +
                   site, data = scaled)
summary(mmdoc.d.temp) # 0.797
AIC(mmdoc.d.temp) # 78.02 tiny bit worse
pred_r_squared(mmdoc.d.temp) #.724 better!

CH4.mod <- mmdoc.d.temp
r.squaredGLMM(CH4.mod)
pred_r_squared(CH4.mod)
summary(CH4.mod)
mmdoc.d.temp <- lm(CH4.ugL ~ logQ + watertemp_C + GPP + ER + DO.obs +
                   site-1, data = scaled)
summary(mmdoc.d.temp)


# CH4.flux ####
mm <- lm(CH4.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER + DO.obs +
           doc_mgl + site, data = scaled)
stepAIC(mm)
Anova(mm)
summary(mm) # .752
AIC(mm) # 91.04

# try DO
mmDO <- lm(CH4.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER + doc_mgl +
            site, data = scaled)
anova(mmDO)
summary(mmDO) # .758
AIC(mmDO) #89.04, better

# try doc
mmDO.doc <- lm(CH4.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER +
                  site, data = scaled)
summary(mmDO.doc) # .762
AIC(mmDO.doc) #87.33, better
pred_r_squared(mmDO.doc) #.676
#try temp
mmDO.doc.temp <- lm(CH4.flux_ugld ~ logQ + depth + GPP + ER + site,
                    data = scaled)
summary(mmDO.doc.temp) # .753
AIC(mmDO.doc.temp) #88.61 worse
pred_r_squared(mmDO.doc.temp) #.666
vif(mmDO.doc) # depth Q and site have high vif - look for correlations
# try depth

mmDO.doc.d <- lm(CH4.flux_ugld ~ logQ + watertemp_C + GPP + ER  + site,
                 data = scaled)
summary(mmDO.doc.d) # .735
pred_r_squared(mmDO.doc.d) # 0.641
AIC(mmDO.doc.d) #92.39 worse, but check vif
vif(mmDO.doc.d) # looks a lot better, still test removing Q or site

mmDO.doc.Q <- lm(CH4.flux_ugld ~ watertemp_C + depth + GPP + ER +
              site, data = scaled)
summary(mmDO.doc.Q) # .715
AIC(mmDO.doc.Q) # 96.32 worse

mmDO.doc.site <- lm(CH4.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER,
                    data = scaled)
summary(mmDO.doc.site) # .62
AIC(mmDO.doc.site) #108.08 much worse

# try gpp
mmDO.doc.gpp <- lm(CH4.flux_ugld ~ logQ + watertemp_C + depth + ER  + site,
                 data = scaled)
summary(mmDO.doc.gpp) # .718
pred_r_squared(mmDO.doc.gpp) # 0.631
AIC(mmDO.doc.gpp) #95.8 worse,
vif(mmDO.doc.gpp) # still bad.

# try gpp with depth removed
mmDO.doc.d.gpp <- lm(CH4.flux_ugld ~ logQ + watertemp_C + ER  + site,
                 data = scaled)
summary(mmDO.doc.d.gpp) # .719
pred_r_squared(mmDO.doc.d.gpp) # 0.636, somewhat better relative to r2, but not really
AIC(mmDO.doc.d.gpp) #94.89 worse,

# depth seems like the best one to take out.
CH4flux.mod <- mmDO.doc.d
r.squaredGLMM(CH4flux.mod)
pred_r_squared(CH4flux.mod)
mmDO.doc.d <- lm(CH4.flux_ugld ~ logQ + watertemp_C + GPP + ER  + site-1,
                 data = scaled)

summary(mmDO.doc.d) # .735

# N2O conc ####
mm <- lm(N2O.ugL ~ logQ + watertemp_C + depth + GPP + ER +
           DO.obs + no3n_mgl + doc_mgl, data = scaled)
Anova(mm)
summary(mm) #.41
AIC(mm) # 131.7 this is after removing site, which none were signifcant

mmgpp <- lm(N2O.ugL ~ logQ + watertemp_C + depth + ER + DO.obs +
              no3n_mgl + doc_mgl, data = scaled)
summary(mmgpp) # .416
AIC(mmgpp) #130.3 a little better
# remove water temp
mmgpp.temp <- lm(N2O.ugL ~ logQ + depth + ER + DO.obs + no3n_mgl + doc_mgl,
              data = scaled)
summary(mmgpp.temp) # .418
AIC(mmgpp.temp) #129.3 a little better
#next DO
mmgpp.temp.DO <- lm(N2O.ugL ~ logQ + depth + ER + no3n_mgl + doc_mgl,
                 data = scaled)
summary(mmgpp.temp.DO) # .426
AIC(mmgpp.temp.DO) #127.7 a little better
pred_r_squared(mmgpp.temp.DO) # 0.32
# next depth
mmgpp.temp.DO.depth <- lm(N2O.ugL ~ logQ + ER + no3n_mgl + doc_mgl,
                          data = scaled)
summary(mmgpp.temp.DO.depth) # .409
AIC(mmgpp.temp.DO.depth) #128.4 slightly worse
# ER
mmgpp.temp.DO.er <- lm(N2O.ugL ~ logQ + depth + no3n_mgl + doc_mgl,
                      data = scaled)
summary(mmgpp.temp.DO.er) # .408
AIC(mmgpp.temp.DO.er) #128.4 worse
# doc
mmgpp.temp.DO.doc <- lm(N2O.ugL ~ logQ + depth + ER + no3n_mgl,
                      data = scaled)
summary(mmgpp.temp.DO.doc) # .407
AIC(mmgpp.temp.DO.doc) #128.5 also worse.

# try combinations: depth and ER:
mmgpp.temp.DO.depther <- lm(N2O.ugL ~ logQ + no3n_mgl + doc_mgl,
                          data = scaled)
summary(mmgpp.temp.DO.depther) # .395
AIC(mmgpp.temp.DO.depther) #128.8 slightly worse
pred_r_squared(mmgpp.temp.DO.depther) # 0.324
mmgpp.temp.DO.depthdoc <- lm(N2O.ugL ~ logQ + no3n_mgl + ER,
                          data = scaled)
summary(mmgpp.temp.DO.depthdoc) # .404
AIC(mmgpp.temp.DO.depthdoc) #127.95 same
pred_r_squared(mmgpp.temp.DO.depthdoc) # 0.342 Better!!
mmgpp.temp.DO.erdoc <- lm(N2O.ugL ~ logQ + no3n_mgl + depth,
                          data = scaled)
summary(mmgpp.temp.DO.erdoc) # .403
AIC(mmgpp.temp.DO.erdoc) #128.03 same
pred_r_squared(mmgpp.temp.DO.erdoc) # 0.351 Better!!
mmgpp.temp.DO.erdocdepth <- lm(N2O.ugL ~ logQ + no3n_mgl,
                          data = scaled)
summary(mmgpp.temp.DO.erdocdepth) # .400
AIC(mmgpp.temp.DO.erdocdepth) #127.44 best
pred_r_squared(mmgpp.temp.DO.erdocdepth) # 0.361 Best!!


vif(mmgpp.temp.DO.erdocdepth) # looks good
N2O.mod <- mmgpp.temp.DO.erdocdepth
r.squaredGLMM(N2O.mod)
pred_r_squared(N2O.mod)

# N2o flux ####
# removed site from beginning because leaving it in
# results in a lower pred r2 in the end
mm <- lm(N2O.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER +
           DO.obs + no3n_mgl + doc_mgl, data = scaled)
Anova(mm)
summary(mm) # .249
AIC(mm) #144.73
# take out ER
mmer <- lm(N2O.flux_ugld ~ logQ + watertemp_C + depth + GPP + DO.obs +
               no3n_mgl + doc_mgl, data = scaled)
summary(mmer) # .266
AIC(mmer) # 142.7, improvement

#take out temp
mmer.temp <- lm(N2O.flux_ugld ~ logQ + depth + GPP + DO.obs +
                   no3n_mgl + doc_mgl, data = scaled)
summary(mmer.temp) # .28
AIC(mmer.temp) # 140.8, improvement

# take out doc
mmer.temp.doc <- lm(N2O.flux_ugld ~ logQ + depth + GPP + DO.obs +
                      no3n_mgl, data = scaled)
summary(mmer.temp.doc) # .295
AIC(mmer.temp.doc) # 138.9, improvement
# take out DO
mmer.temp.doc.DO <- lm(N2O.flux_ugld ~ logQ + depth + GPP +
                      no3n_mgl, data = scaled)
summary(mmer.temp.doc.DO) # .300
AIC(mmer.temp.doc.DO) # 137.714, improvement
# take out Q
mmer.temp.doc.DO.Q <- lm(N2O.flux_ugld ~ depth + GPP + no3n_mgl, data = scaled)
summary(mmer.temp.doc.DO.Q) # .287
AIC(mmer.temp.doc.DO.Q) # 137.72, almost the same

pred_r_squared(mmer.temp.doc.DO.Q)
# GPP
mmer.temp.doc.DO.Q.gpp <- lm(N2O.flux_ugld ~ depth + no3n_mgl, data = scaled)
summary(mmer.temp.doc.DO.Q.gpp) # .289
AIC(mmer.temp.doc.DO.Q.gpp) # 136.7, better

vif(mmer.temp.doc.DO.Q.gpp)
pred_r_squared(mmer.temp.doc.DO.Q.gpp)

N2Oflux.mod <- mmer.temp.doc.DO.Q.gpp
r.squaredGLMM(N2Oflux.mod)
pred_r_squared(N2Oflux.mod)

# O2 conc ####
mm <- lm(DO.obs ~ logQ + watertemp_C + depth + GPP + ER + doc_mgl + no3n_mgl +
           site, data = scaled)
summary(mm) #.715
AIC(mm) #95.99
Anova(mm)
# remove discharge
mmQ <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + doc_mgl +
               no3n_mgl + site, data = scaled)
summary(mmQ) #.722
AIC(mmQ) #94.01 better
Anova(mmQ)
# remove no3
mmQ.no3 <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + doc_mgl +
           site, data = scaled)
summary(mmQ.no3) #.718
AIC(mmQ.no3) #96.4 worse
# remove doc
mmQ.doc <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + no3n_mgl +
           site, data = scaled)
summary(mmQ.doc) #.705
AIC(mmQ.doc) #96.4 worse
vif(mmQ)
# remove site
mmQ.site <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + doc_mgl + no3n_mgl,
          data = scaled)
AIC(mmQ.site) #94.8
summary(mmQ.site)# .696
vif(mmQ.site)
# remove depth
mmQ.d <- lm(DO.obs ~ watertemp_C  + GPP + ER + doc_mgl + no3n_mgl + site,
          data = scaled)
AIC(mmQ.d) #97.8
summary(mmQ.d)# .696 # worse
vif(mmQ.site)
# try removing combinations of drivers:
# site and no3
mmQ.siteno3 <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + doc_mgl,
                  data = scaled)
summary(mmQ.siteno3) #.699
AIC(mmQ.siteno3) #96.0
pred_r_squared(mmQ.siteno3) #.670

mmQ.no3doc <- lm(DO.obs ~ watertemp_C + depth + GPP + ER +
             site, data = scaled)
summary(mmQ.no3doc) #.713
AIC(mmQ.no3doc) #96.78
pred_r_squared(mmQ.no3doc) #.657
mmQ.sitedoc <- lm(DO.obs ~ watertemp_C + depth + GPP + ER + no3n_mgl,
                  data = scaled)
summary(mmQ.sitedoc) #.703
AIC(mmQ.sitedoc) #92.8 better
pred_r_squared(mmQ.sitedoc) #.667
Anova(mmQ.sitedoc)

#remove no3
mmQ.sitedocno3 <- lm(DO.obs ~ watertemp_C + depth + GPP + ER,
                  data = scaled)
summary(mmQ.sitedocno3) #.705 same
AIC(mmQ.sitedocno3) #94.01 little worse
pred_r_squared(mmQ.sitedocno3) #.676 little better.
# remove depth
mmQ.sitedoc.depth <- lm(DO.obs ~ watertemp_C + GPP + ER + no3n_mgl,
                  data = scaled)
summary(mmQ.sitedoc.depth) #.700
AIC(mmQ.sitedoc.depth) #92.3 little better
pred_r_squared(mmQ.sitedoc.depth) #.669 little better

#remove both
mmQ.sitedoc.depthno3 <- lm(DO.obs ~ watertemp_C + GPP + ER,
                  data = scaled)
summary(mmQ.sitedoc.depthno3) #.701
AIC(mmQ.sitedoc.depthno3) #93.86 little better
pred_r_squared(mmQ.sitedoc.depthno3) #.674 little better
Anova(mmQ.sitedoc.depthno3)

O2.mod <- mmQ.sitedoc.depthno3
r.squaredGLMM(O2.mod)
pred_r_squared(O2.mod)
summary(O2.mod)

# O2 flux ####
mm <- lm(O2.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER + doc_mgl +
           no3n_mgl +site, data = scaled)
summary(mm) #.911
stepAIC(mm) #35.0
Anova(mm)
# remove site
mmsite <- lm(O2.flux_ugld ~ logQ + watertemp_C + depth + GPP + ER + doc_mgl +
               no3n_mgl, data = scaled)
summary(mmsite) #.918
AIC(mmsite) #26.98 big improvement
# remove Q
mmsite.Q <- lm(O2.flux_ugld ~ watertemp_C + depth + GPP + ER + doc_mgl + no3n_mgl,
               data = scaled)
summary(mmsite.Q) #.919
AIC(mmsite.Q) #25.2 better
# remove doc
mmsite.Q.doc <- lm(O2.flux_ugld ~ watertemp_C + depth + GPP + ER + no3n_mgl,
                   data = scaled)
summary(mmsite.Q.doc) #.922
AIC(mmsite.Q.doc) #23.2 better
# remove temp
mmsite.Q.doc.temp <- lm(O2.flux_ugld ~ depth + GPP + ER + no3n_mgl, data = scaled)
summary(mmsite.Q.doc.temp) #.922
AIC(mmsite.Q.doc.temp) #22.03 better
Anova(mmsite.Q.doc.temp)
vif(mmsite.Q.doc.temp)

O2flux.mod <- mmsite.Q.doc.temp
r.squaredGLMM(O2flux.mod)
pred_r_squared(O2flux.mod) #0.904


# Compile hand modeled results ####
rsqs <- data.frame()
mod_fits <- data.frame()

handmods <- list(O2.mod, CO2.mod, CH4.mod, N2O.mod,
                 O2flux.mod, CO2flux.mod, CH4flux.mod, N2Oflux.mod)

gasses <- c('N2O', 'O2', 'CO2', 'CH4')
for(i in 1:8){

  if(i <= 4){ flux = FALSE} else { flux = TRUE}
  gas = gasses[(i %% 4) + 1]
  mm <- handmods[[i]]
  out <- get_mod_resultss(mod_fits, out$mod)
}

saveRDS(list(mod_fits = mod_fits,
             rsqs = rsqs), (mm, gas, flux)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_row
        'data/model_fits_steps_rsqs_lm_manual_AIC_elimination.rds')


# Automatic AIC selection ####
rsqs <- data.frame()
mod_fits <- data.frame()
mod_steps <- data.frame()

# CO2
mm <- lm(CO2.ugL ~ logQ + watertemp_C + slope_nhd + depth + GPP + ER + DO.obs +
           doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = "CO2")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lm(CO2.flux_ugld ~ logQ + watertemp_C + slope_nhd + depth + GPP + ER +
           DO.obs + doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = 'CO2', flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

# CH4
mm <- lm(CH4.ugL ~ logQ + watertemp_C + slope_nhd + depth+ GPP + ER + DO.obs +
           doc_mgl + CO2.ugL + site, data = scaled)
  out <- get_mod_results(mm, gas = 'CH4')
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

mm <- lm(CH4.flux_ugld ~ logQ + watertemp_C + slope_nhd  + GPP + ER + DO.obs +
           doc_mgl + CO2.ugL + site, data = scaled)
  out <- get_mod_results(mm, gas = 'CH4', flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

# N2O
mm <- lm(N2O.ugL ~ logQ + watertemp_C + slope_nhd + depth +GPP + ER +
           DO.obs + no3n_mgl + doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = 'N2O')
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

mm <- lm(N2O.flux_ugld ~ logQ + watertemp_C + slope_nhd + depth + GPP + ER +
           DO.obs + no3n_mgl + doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = "N2O", flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

# O2
mm <- lm(DO.obs ~ logQ + watertemp_C + slope_nhd + depth +GPP + ER +
           doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = "O2")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lm(O2.flux_ugld ~ logQ + watertemp_C + slope_nhd + depth + GPP + ER +
           doc_mgl + site, data = scaled)
  out <- get_mod_results(mm, gas = "O2", flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

saveRDS(list(mod_fits = mod_fits,
             mod_steps = mod_steps,
             rsqs = rsqs),
        'data/model_fits_steps_rsqs_lm_AIC_elimination.rds')


# Analyze water chem data ####
spchem <- read_csv("data/water_chemistry/all_grab_data.csv") %>%
  filter(siteID %in% c("NHC", "UNHC")) %>%
  dplyr::select(-flagID, -flagComment, -methodDetail, -writeInMethod, -regionID, -method) %>%
  group_by(siteID, dateTimeUTC, variable) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%#data.frame()
  filter(!is.na(siteID), !is.na(dateTimeUTC),
         !(variable %in% c('Br'))) %>%
  rename(site = siteID, DateTime_UTC = dateTimeUTC)

# pair with discharge and temperature
dat <- read_csv('../NHC_2019_metabolism/data/rating_curves/interpolatedQ_allsites_modified.csv') %>%
  rename(NHC = NHC.Q, UNHC = UNHC.Q) %>%
  dplyr::select(DateTime_UTC, NHC, UNHC) %>%
  pivot_longer(cols = -DateTime_UTC, names_to = "site", values_to = 'discharge')

chem <- spchem %>%
  left_join(dat, by = c('site', 'DateTime_UTC')) %>%
  mutate(doy = as.numeric(format(DateTime_UTC, '%j')))
ggplot(chem, aes(log(discharge), value, col = site)) +
  geom_point()+
  geom_smooth(method = lm) +
  facet_wrap(.~variable, scales = "free_y")

# chem.pca <- chem %>%
#   select(-siteID, -dateTimeUTC, -discharge, -TOC, -temp.water) %>%
#   prcomp()
#
# png("figures/gas/waterchem_pca.png", height = 5, width = 5, units = "in", res = 300)
# autoplot(chem.pca, data = chem, colour = "siteID", size = 2,
#          loadings = TRUE, loadings.label.colour = 'black',
#          loadings.label = TRUE, loadings.label.size = 5) +
#   # scale_color_gradient(low = "black", high = "red") +
#   labs(title = "Water Chem by Discharge @ NHC, UNHC")
# dev.off()
# biplot(chem.pca, pch = 20)
