# GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(MASS)
library(HH)
library(MuMIn)
# library(lme4)
library(nlme)
library(lmerTest)
library(car)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe_individual_samples.csv")

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
         log_no3n = log(no3n_mgl),
         log_doc = log(doc_mgl)
         ) %>%
    dplyr::select(site, sample, logQ, watertemp_C, ER, GPP, DO.obs, slope_mm, depth,
            log_no3n, log_doc, ends_with(c("ugld", "ugL", 'mgm2d'))) %>%
    mutate(across(-any_of(c('site', 'sample')), ~ scale(.)[,1, drop = T]),
         across(all_of(c('site', 'sample')), ~ factor(.)))
# scaled <- dat %>%
#     filter(site != 'MC751',
#          !is.na(CH4.ugL),
#          !is.na(GPP)) %>%
#     mutate(logQ = log(discharge),
#          NER = ER - GPP,
#          DO.persat = DO.obs/DO.sat,
#          # logWRT = log(1000*depth*width_march_m/discharge),
#          no3n_mgl = ifelse(no3n_mgl == 0, 0.0015, no3n_mgl), # replace zero with mdl
#          site = factor(site, levels=c('UNHC','WBP','WB','CBP','PM','NHC')),
#          log_no3n = log(no3n_mgl),
#          log_doc = log(doc_mgl),
#          N2O.ugL = case_when(N2O.ugL == 0 ~ 0.015,
#                              TRUE ~ N2O.ugL)) %>%
#     mutate(across((ends_with('ugL')), ~log(.x)))%>%
#     dplyr::select(site, sample, logQ, watertemp_C, ER, GPP, DO.obs, slope_mm, depth,
#             log_no3n, log_doc, ends_with(c("ugld", "ugL", 'mgm2d'))) %>%
#     mutate(across(-any_of(c('site', 'sample')), ~ scale(.)[,1, drop = T]),
#          across(all_of(c('site', 'sample')), ~ factor(.)))
#

preds <- scaled %>%
  dplyr::select(site, sample, logQ, watertemp_C, GPP, ER, slope_mm, DO.obs,
         log_no3n, log_doc, depth)

pred_cov <- data.frame(cov(preds[,3:11]))
write_csv(pred_cov, 'data/linear_models/predictor_covariance_matrix.csv')

# correlated predictors(r > 0.5)
cor.test(preds$logQ, preds$GPP, method = 'pearson')
cor.test(preds$logQ, preds$log_doc, method = 'pearson')
cor.test(preds$logQ, preds$depth, method = 'pearson')
cor.test(preds$watertemp_C, preds$GPP, method = 'pearson')
cor.test(preds$ER, preds$DO.obs, method = 'pearson')

# functions for assessing models ####
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
determine_site_signif <- function(y, scaled, flux = FALSE){
    scaled$y <- y
    # how much additional variation is explained by site, once slope and depth are acct'd for?
    grp <- scaled %>%
        group_by(sample, site) %>%
        summarize(across(any_of(c('y', 'slope_mm', 'depth')), mean, na.rm = T))
    lme0site <- lme4::lmer(y~ (1|site) + slope_mm + depth, data=grp)
    if(flux){lme0site <- lme4::lmer(y ~ (1|site) + slope_mm, data=grp)}
    vnc <- as.data.frame(VarCorr(lme0site))$vcov
    add_var <- vnc[1]/sum(vnc)


    av <- anova(lm(y ~ site, data = grp))

    p_val <- av$`Pr(>F)`[1]

    return(data.frame(add_var = add_var,
                      p_val = p_val))
}

# note, a singular fit just means one of your variances is very close to zero, which we expect in some of these models!
site_sig <- data.frame()
site_sig <- determine_site_signif(scaled$CO2.ugL, scaled) %>%
    mutate(gas = 'CO2conc') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$CO2.flux_mgm2d, scaled, TRUE) %>%
    mutate(gas = 'CO2flux') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$DO.obs, scaled) %>%
    mutate(gas = 'O2conc') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$O2.flux_mgm2d, scaled, TRUE) %>%
    mutate(gas = 'O2flux') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$CH4.ugL, scaled) %>%
    mutate(gas = 'CH4conc') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$CH4.flux_mgm2d, scaled, TRUE) %>%
    mutate(gas = 'CH4flux') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$N2O.ugL, scaled) %>%
    mutate(gas = 'N2Oconc') %>%
    bind_rows(site_sig)
site_sig <- determine_site_signif(scaled$N2O.flux_mgm2d, scaled, TRUE) %>%
    mutate(gas = 'N2Oflux') %>%
    bind_rows(site_sig)

calc_loocv_rmse <- function(scaled, formula){
    rsqe <- rep(NA, nrow(scaled))
    for(i in 1:nrow(scaled)){
        sc <- scaled[-i,]
        mm <- lmer(formula, data = sc)
        p <- try( predict(mm, newdata = scaled, allow.new.levels = TRUE)[i] )
        if(inherits(p, 'try-error')) next

        rsqe[i] <- (p - scaled$CH4.ugL[i])^2
    }
    loocv = sqrt(mean(rsqe, na.rm = T))
    return(loocv)
}
search_lmer <- function(y, preds, gas, flux = 'no'){
    if(gas == 'O2'){
        preds <- dplyr::select(preds, -DO.obs)
    }

    vars <- colnames(preds)[3:ncol(preds)]
    nvar <- min(length(vars), 5)
    preds$y <- y
    m1 <- lmer(y ~ (1|sample) + (1|site), data = preds)
    mods <- data.frame(sample = TRUE, site = TRUE) %>%
        mutate(aicc = AICc(m1),
               singular = isSingular(m1),
               vif = 0)
    mods <- bind_cols(mods, r.squaredGLMM(m1))
    # models without site:
    for(nv in 1:nvar){
        vsets <- combn(length(vars), nv)
        for(i in 1:ncol(vsets)){
            vv = vars[vsets[,i]]
            r <- as.data.frame(matrix(ncol = nv, rep(TRUE,nv)))
            colnames(r) <- vv
            m1 <- lmer(paste0('y ~ (1|sample) + ', paste(vv, collapse = ' + ')),
                       data = preds)
            r1 <- mutate(r, sample = TRUE)
            m2 <- lmer(paste0('y ~ (1|sample) + (1|site) +',
                              paste(vv, collapse = ' + ')),
                       data = preds)
            r2 <- mutate(r1, site = TRUE)
            r1$aicc <- AICc(m1)
            r1$singular = isSingular(m1)
            r1 <- bind_cols(r1, r.squaredGLMM(m1))
            r1$rmse <- sqrt(mean((residuals(m1))^2))
            r2$aicc <- AICc(m2)
            r2 <- bind_cols(r2, r.squaredGLMM(m2))
            r2$singular = isSingular(m2)
            r2$rmse <- sqrt(mean((residuals(m2))^2))
            r1$vif <- r2$vif <- 0
            if(nv >1){
                r1$vif <- max(vif(m1))
                r2$vif <- max(vif(m2))
            }
            mods <- bind_rows(mods, r2, r1)
        }
    }

    if(flux == 'flux'){
        mods <- filter(mods, is.na(depth)) %>%
            dplyr::select(-depth)
    }

    if(gas != 'N2O'){
        mods <- filter(mods, is.na(log_no3n)) %>%
            dplyr::select(-log_no3n)
    }

    mods <- mods %>% tibble() %>%
        filter(vif < 5,
               !(!is.na(log_doc) & !is.na(logQ))) %>%
        arrange(aicc) %>%
        mutate(delta_aicc = aicc - min(aicc)) %>%
        # filter(delta_aicc < 5) %>%
        slice(1:10) %>%
        mutate(rel_likelihood = exp(-0.5 * delta_aicc),
               aicc_weight = rel_likelihood/sum(rel_likelihood),
               gas = gas,
               flux = flux)

    site = TRUE
    if(is.na(mods[1,'site'])) site = FALSE

    m <- mods[1,] %>%
        dplyr::select(-site, -sample) %>%
        select_if(isTRUE)
    vars <- colnames(m)

    mm <- lmer(paste0('y ~ (1|sample) + ', paste(vars, collapse = ' + ')),
               data = preds)
    if(site){
        mm <- lmer(paste0('y ~ (1|sample) + (1|site) + ',
                          paste(vars, collapse = ' + ')),
               data = preds)
    }

    return(list(mods = mods, m1 = mm))
}

best_lmes <- data.frame()
out_ch4 <- search_lmer(scaled$CH4.ugL, preds, 'CH4')
best_lmes <- bind_rows(best_lmes, out_ch4$mods)
out_ch4f <- search_lmer(scaled$CH4.flux_mgm2d, preds, 'CH4', 'flux')
best_lmes <- bind_rows(best_lmes, out_ch4f$mods)
out_N2O <- search_lmer(scaled$N2O.ugL, preds, 'N2O')
best_lmes <- bind_rows(best_lmes, out_N2O$mods)
out_N2Of <- search_lmer(scaled$N2O.flux_mgm2d, preds, 'N2O', 'flux')
best_lmes <- bind_rows(best_lmes, out_N2Of$mods)
out_CO2 <- search_lmer(scaled$CO2.ugL, preds, 'CO2')
best_lmes <- bind_rows(best_lmes, out_CO2$mods)
out_CO2f <- search_lmer(scaled$CO2.flux_mgm2d, preds, 'CO2', 'flux')
best_lmes <- bind_rows(best_lmes, out_CO2f$mods)

# out_O2 <- search_lmer(scaled$DO.obs, preds, 'O2')
# best_lmes <- bind_rows(best_lmes, out_O2$mods)
# out_O2f <- search_lmer(scaled$O2.flux_mgm2d, preds, 'O2', 'flux')
# best_lmes <- bind_rows(best_lmes, out_O2f$mods)
write_csv(best_lmes, 'data/linear_models/best_lme_summaries.csv')

mods <- list(CH4.conc = out_ch4$m1,
             CH4.flux = out_ch4f$m1,
             CO2.conc = out_CO2$m1,
             CO2.flux = out_CO2f$m1,
             N2O.conc = out_N2O$m1,
             N2O.flux = out_N2Of$m1)

saveRDS(mods, 'data/linear_models/best_lmes.rds')

summary(out_CO2f$m1)

#one thought was that GHG patterns by site were more stable than between sites
#so should probably test the site as a fixed effect
grouped <- scaled1 %>% group_by(site, sample) %>%
    summarize(across(ends_with(c('ugL', 'mgm2d')), ~mean(.x, na.rm = T)))
summary(aov(CO2.ugL ~ site, data = grouped))         # not significant
summary(aov(CO2.flux_mgm2d ~ site, data = grouped))   # not significant
summary(aov(CH4.ugL ~ site, data = grouped))         # *** p = 0.000758
summary(aov(CH4.flux_mgm2d ~ site, data = grouped))   # not significant
summary(aov(N2O.ugL ~ site, data = grouped))         # not significant
summary(aov(N2O.flux_mgm2d ~ site, data = grouped))   # not significant

summary(lm(CH4.ugL ~ site-1, data = grouped))

# this method is no longer in use 6/2022
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
