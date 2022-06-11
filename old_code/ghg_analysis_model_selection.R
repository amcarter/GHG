# GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(lme4)
# library(lmerTest)
library(glmulti)
library(leaps)

library(MASS)
library(HH)
library(MuMIn)

setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")

dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
# d2 <- read_csv("data/ghg_flux_complete_drivers_dataframe_noNAs.csv")
filter(dat, site !='MC751') %>% 
  mutate(date = as.Date(group)) %>%
ggplot(aes(date, doc_mgl, col = site)) +
  geom_line() + 
  geom_point() +
  facet_wrap(~site)

# linear mixed effects models ####
# rescale covariates, normalize each to the mean
scaled <- dat %>%
  filter(site != 'MC751',
         !is.na(CH4.ugL),
         !is.na(GPP)) %>%
  mutate(logQ = log(discharge),
         NER = ER - GPP,
         DO.persat = DO.obs/DO.sat,
         log_no3 = case_when(no3n_mgl== 0 ~ log(0.0015),
                                   TRUE ~ log(no3n_mgl))) %>% 
         # no3n_mgl = ifelse(no3n_mgl > 0.6, NA, no3n_mgl)) %>%# remove high NO3, it is too high leverage
  dplyr::select(site, habitat, logQ, watertemp_C, ER, GPP, NER, K600,
         DO.obs, DO.persat, slope_wbx,
         log_no3, nh4n_mgl, doc_mgl, ends_with("ugld"), ends_with("ugL")) %>%
  mutate(across(-c(site, habitat), ~ scale(.)[,1, drop = T]), 
         across(c(site, habitat), ~ factor(.)))

par(mfrow = c(3,4))
for(i in 3:21){
  
  plot(density(scaled[,i, drop = T], na.rm = T))}
preds <- scaled %>%
  dplyr::select(logQ, watertemp_C, GPP, ER, K600, slope_wbx, DO.obs, 
         log_no3, CO2.ugL, doc_mgl) 

cov(preds)
# There is a lot of covariance between:
# O2, and ER. I will let ER represent both of them
# CO2 and GPP and ER. I will exclude CO2
# DOC and logQ. I wll exclude DOC

# run an exhaustive set of LMEs using the Leaps Package
get_mod_results <- function(mm, gas, flux = F){
  m <- stepAIC(mm, direction = 'both', trace = F)
  print("Variance inflation factors - should be less than ~5")
  print(vif(m))
  r <- data.frame(r.squaredGLMM(m)) %>%
    mutate(gas = gas, flux = flux,
           formula = as.character(summary(m)$call)[2])
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
  
  steps <- m$anova %>%
    as_tibble() %>%
    rename(Resid_Df = 'Resid. Df', resid_dev = 'Resid. Dev') %>%
    mutate(effect = "fixed",
           gas = gas,
           flux = flux) 
  
  
  return(list(r2 = r,
              steps = steps,
              mod = ci))
}

rsqs <- data.frame()
mod_fits <- data.frame()
mod_steps <- data.frame()

# CO2####
# test if site should be included as a random effectCO2
nulm <- lm(CO2.flux_ugld ~ slope_wbx, data = scaled)
nm <- lm(CO2.flux_ugld ~ slope_wbx + K600, data = scaled)
nmb <- lm(CO2.flux_ugld ~ K600, data = scaled)
summary(nmb) 
AICc(mix)
mix <- lme4::lmer(CO2.flux_ugld ~ (1|site), data = scaled)
summary(mix)  # 0 % additional variance explained

# likelihood ratio test: I'm just not sure how to interpret the output
# my_stat = 2*(logLik(mix) - logLik(nulm,REML = TRUE))
# LRT_stat = numeric(1000)
# for(i in 1:1000){
#   y = simulate(nulm)$sim_1
#   null_md = lm(y ~ slope_wbx, data = scaled)
#   mixed_md = lmer(y ~ slope_wbx + (1|site), data = scaled)
#   LRT_stat[i] = as.numeric(2*(logLik(mixed_md) - logLik(null_md,REML = TRUE)))
# }
# sum(LRT_stat > my_stat)/1000

mm <- lm(CO2.ugL ~ slope_wbx, data = scaled)
out <- get_mod_results(mm, gas = "CO2")
rsqs <- bind_rows(rsqs, out$r2)
mod_steps <- bind_rows(mod_steps, out$steps)
mod_fits <- bind_rows(mod_fits, out$mod)
Q + watertemp_C + slope_nhd + GPP + ER + DO.obs +
           doc_mgl, data = scaled)
mm <- lm(CO2.flux_ugld ~ logQ + watertemp_C + slope_nhd + GPP + ER + 
           DO.obs + doc_mgl, data = scaled)
out <- get_mod_results(mm, gas = 'CO2', flux = TRUE)
rsqs <- bind_rows(rsqs, out$r2)
mod_steps <- bind_rows(mod_steps, out$steps)
mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lme4::lmer(CO2.flux_ugld ~ slope_wbx + (1|site), data = scaled)
summary(mm)
0.05/(0.05 + 0.985) # 4.8% additional variance explained

mm <- lme4::lmer(CH4.ugL ~ slope_wbx + (1|site), data = scaled)
summary(mm)  # 0 % additional variance explained

mm <- lme4::lmer(CH4.flux_ugld ~ slope_wbx + (1|site), data = scaled)
summary(mm)
0.123/(0.123 + 0.802) # 13% additional variance explained

mm <- lme4::lmer(N2O.ugL ~ slope_wbx + (1|site), data = scaled)
summary(mm)  # 0 % additional variance explained

mm <- lme4::lmer(N2O.flux_ugld ~ slope_wbx + (1|site), data = scaled)
summary(mm)
0.036/(0.036 + 0.828) # 4.1% additional variance explained

# it doesn't seem worth it to include site as a random effect in my models.
# Use leaps package to test all compination of plain lms

  
# N2O
mm <- lm(N2O.ugL ~ logQ + watertemp_C + slope_wbx + GPP + ER +
           DO.obs + no3n_mgl + doc_mgl, data = scaled)
  out <- get_mod_results(mm, gas = 'N2O')
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

mm <- lm(N2O.flux_ugld ~ logQ + watertemp_C + slope_wbx + GPP + ER +
           DO.obs + no3n_mgl + doc_mgl, data = scaled)
  out <- get_mod_results(mm, gas = "N2O", flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)
  
# O2
mm <- lm(DO.obs ~ logQ + watertemp_C + slope_wbx + GPP + ER + 
           doc_mgl, data = scaled)
  out <- get_mod_results(mm, gas = "O2")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lm(O2.flux_ugld ~ logQ + watertemp_C + slope_wbx + GPP + ER + 
           doc_mgl, data = scaled)
  out <- get_mod_results(mm, gas = "O2", flux = TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_steps <- bind_rows(mod_steps, out$steps)
  mod_fits <- bind_rows(mod_fits, out$mod)

saveRDS(list(mod_fits = mod_fits,
             mod_steps = mod_steps,
             rsqs = rsqs), 
        'data/model_fits_steps_rsqs_lm_AIC_elimination.rds')

# # LMER versions ####  
# I think this code won't work, even if you uncomment it. 
# In general though, my LMER approach was a model like this:
 mm <- lmer(CO2.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs + doc_mgl +
               (1|site), data = scaled)
# Then using step from the lmerTest package


get_mod_results_rand <- function(mm, gas, flux = F){
  mm_step <- lmerTest::step(mm)
  m <- get_model(mm_step)
  r <- data.frame(r.squaredGLMM(m)) %>%
    mutate(gas = gas, flux = flux,
           formula = as.character(summary(m)$call)[2])
  cf <- summary(m)$coefficients %>%
    as_tibble() %>%
    mutate(predictor = rownames(summary(m)$coefficients))%>%
    rename(estimate = Estimate, std_error = 'Std. Error', t_val = 't value',
           pr_greaterthan_t = 'Pr(>|t|)')
  ci <- data.frame(confint(m)) %>%
    mutate(predictor = rownames(confint(m)))
  steps_r <- mm_step$random %>%
    as_tibble() %>%
    mutate(predictor = rownames(mm_step$random),
           effect = "random") %>%
    rename(pr_greaterthan_chisq = 'Pr(>Chisq)')
  steps <- mm_step$fixed %>%
    as_tibble() %>%
    mutate(predictor = rownames(mm_step$fixed),
           effect = "fixed") %>%
    rename_all(dplyr::recode, 'Sum of Sq'='Sum Sq') %>%
    rename(sum_sq = 'Sum Sq',  F_val = 'F value',
           pr_greaterthan_F = 'Pr(>F)') %>%
    full_join(cf, by = 'predictor') %>%
    left_join(ci, by = 'predictor') %>%
    bind_rows(steps_r) %>%
    mutate(gas = gas,
           flux = flux)

  return(list(r2 = r,
              mod = steps))
}

print("Variance inflation factors - should be less than ~5")
  v <-vif(mm)
  vf <- tibble(predictor = names(v), vif = v) %>%
          mutate(gas = gas, flux = flux )

function(formula, dat, gas){
  mm <- lme4::lmer(formula, data = dat)
  r <- data.frame(r.squaredGLMM(mm)) %>%
    mutate(gas = gas, flux = flux,
           form = as.character(formula))
  # calculate predictive r2
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
  
  steps <- m$anova %>%
    as_tibble() %>%
    rename(Resid_Df = 'Resid. Df', resid_dev = 'Resid. Dev') %>%
    mutate(effect = "fixed",
           gas = gas,
           flux = flux) 
   
  
}

mod_fits <- data.frame()
rsqs <- data.frame()

mm <- lme4::lmer(CO2.ugL ~ logQ + watertemp_C + GPP + ER +
             slope_wbx + (1|site), data = scaled)
  out <- get_mod_results_rand(mm, gas = "CO2")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lmer(CO2.flux_ugld ~ logQ + watertemp_C + GPP + ER +
             slope_wbx + (1|site), data = scaled)
  out <- get_mod_results_rand(mm, gas = "CO2", TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)

# CH4
mm <- lmer(CH4.ugL ~ logQ + watertemp_C + GPP + ER +
             slope_wbx + (1|site), data = scaled)
  out <- get_mod_results_rand(mm, gas = "CH4")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lmer(CH4.flux_ugld ~ logQ + watertemp_C + GPP + ER + 
             slope_wbx + (1|site), data = scaled)
  out <- get_mod_results_rand(mm, gas = "CH4", TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)

# N2O
scaled_N <- scaled %>% filter(!is.na(no3n_mgl))
mm <- lmer(N2O.ugL ~ logQ + watertemp_C + GPP + ER +
            no3n_mgl + slope_wbx + (1|site), data = scaled_N)
  out <- get_mod_results_rand(mm, gas = "N2O")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)
mm <- lmer(N2O.flux_ugld ~ logQ + watertemp_C + GPP + ER +
            no3n_mgl + slope_wbx + (1|site), data = scaled_N)
  out <- get_mod_results_rand(mm, gas = "N2O", TRUE)
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)

# O2
mm <- lm(DO.obs ~ logQ + watertemp_C + GPP + ER +  doc_mgl + no3n_mgl +
            CO2.ugL, data = scaled)
  out <- get_mod_results(mm, gas = "O2")
  rsqs <- bind_rows(rsqs, out$r2)
  mod_fits <- bind_rows(mod_fits, out$mod)

# write_csv(mod_fits, "data/model_fits_lme_backward_elimination.csv")
# write_csv(rsqs, 'data/model_fits_formula_rsq.csv')
#
# saveRDS(list(mod_fits = mod_fits,
#              mod_steps = mod_steps,
#              rsqs = rsqs),
#         'data/model_fits_steps_rsqs_lm_AIC_elimination_vel_not_Q.rds')

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