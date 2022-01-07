# plot individual regressions based on lme model drivers

library(tidyverse)
library(ggpubr)

mdat <- readRDS('data/model_fits_steps_rsqs_lm_AIC_elimination.rds')
mdat <- readRDS('data/linear_models/model_fits_steps_rsqs_lm_manual_AIC_backward_elimination_03_05.rds')
rsqs <- mdat$rsqs
dvs <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
dvs<- read_csv("data/sites_slope_comparison.csv")%>%
  select(site, slope_deg = slope) %>%
  right_join(dvs, by = 'site')
fits <- mdat$mod_fits %>%
  as_tibble() %>%
  filter(predictor != '(Intercept)') %>%
  rename(CI_2.5 = X2.5.., CI_97.5 = X97.5.., p_val = pr_greaterthan_t)

# write_csv(fits, 'data/model_fits_lm_AIC_elimination.csv')

# lmer
# fits <- read_csv("data/model_fits_lme_backward_elimination.csv")
# rsqs <- read_csv("data/model_fits_formula_rsq.csv")
# fits <- fits %>%
#   # filter(Eliminated == 0, !(predictor %in% c('(Intercept)', '(1 | site)'))) %>%
#   mutate(p_val = case_when(!is.na(pr_greaterthan_t) ~ pr_greaterthan_t,
#                            !is.na(pr_greaterthan_chisq) ~ pr_greaterthan_chisq,
#                            !is.na(pr_greaterthan_F) ~ pr_greaterthan_F)) %>%
#   select(gas, flux, predictor, effect, p_val, estimate, std_error, sum_sq,
#          CI_2.5 = X2.5.., CI_97.5 = X97.5..) 
# 
# write_csv(fits, 'data/model_fits_lme_backward_elimination_condensed.csv')
dat <- dvs %>% 
  filter(!is.na(CH4.ugL),
         !is.na(GPP),
         site != 'MC751')%>%
  mutate(#CO2.ugL2 = CO2.ugL, 
         gas = factor(gas, levels = c('CO2', 'O2', 'CH4', 'N2O')),
         O2.ugL = DO.obs*1000,
         logQ = log(discharge),
         vel_ms = discharge/depth/width_march_m,
         no3n_mgl = ifelse(no3n_mgl >= 0.6, NA, no3n_mgl)) %>%
  pivot_longer(cols = any_of(unique(fits$predictor)), values_to = 'value',
               names_to = 'predictor') %>%
  # rename(CO2.ugL = CO2.ugL2) %>%
  dplyr::select(site, date, group, ends_with(c('.ugL', 'ugld')), 
                predictor, value)

 lm_nums <- function(df){
      m <- lm(gas ~ value, df);
      eq <- paste0("slope = ", format(coef(m)[2], digits = 2), 
                    ', r2 = ', format(summary(m)$r.squared, digits = 3))

      return(eq)
    }

# concentrations
 conc <- fits %>%
  filter(flux == F) %>%
  dplyr::select(-flux, -effect)
 concF <- fits %>%
  filter(flux == T) %>%
  dplyr::select(-flux, -effect)

preds <- dat %>%
  pivot_longer(cols = ends_with(c('.ugL', 'ugld')), names_to = c('gas', 'rate'),
                                values_to = 'gas_value', 
                                names_pattern = '([A-Z0-9]+).([a-zA-Z_]+)') %>%
  mutate(value = case_when(gas == 'O2' & predictor %in%
                       unique(conc$predictor[conc$gas == "O2"]) ~ value,
                     gas == 'O2' & rate == 'flux_ugld' & predictor %in% 
                       unique(concF$predictor[concF$gas == 'O2']) ~ value,
                     gas == 'CO2' & rate == 'ugL' & predictor %in% 
                       unique(conc$predictor[conc$gas == 'CO2']) ~ value,
                     gas == 'CO2' & rate == 'flux_ugld' & predictor %in% 
                       unique(concF$predictor[concF$gas == 'CO2']) ~ value,
                     gas == 'CH4' & rate == 'ugL' & predictor %in% 
                       unique(conc$predictor[conc$gas == 'CH4']) ~ value,
                     gas == 'CH4' & rate == 'flux_ugld' & predictor %in% 
                       unique(concF$predictor[concF$gas == 'CH4']) ~ value,
                     gas == 'N2O' & rate == 'ugL' & predictor %in% 
                       unique(conc$predictor[conc$gas == 'N2O']) ~ value,
                     gas == 'N2O' & rate == 'flux_ugld' & predictor %in% 
                       unique(concF$predictor[concF$gas == 'N2O']) ~ value,
                     TRUE ~ NA_real_),
         rate = ifelse(rate == 'ugL', 'conc', 'flux'),
         gas = factor(gas, levels = c('CO2', 'O2', 'CH4', 'N2O'))) #%>%
  # filter(!(predictor %in% c('doc_mgl', 'DO.obs')))

phys_preds <- filter(preds, predictor %in% 
                       c('watertemp_C', 'slope_deg', 'logQ', 'vel_ms', 'depth'))
bgc_preds <- filter(preds, predictor %in% 
                      c('GPP', 'ER', 'CO2.ugL', 'doc_mgl', 'DO.obs','no3n_mgl'))# %>%
  # mutate(preds = factor(predictor, 
  #                       levels = c('GPP', 'ER', 'CO2.ugL', 'no3n_mgl')))
bgc <- ggplot(bgc_preds, aes(value, gas_value, col = rate)) +
   geom_point() +
   geom_smooth(method = lm)  +
   facet_grid(gas~predictor, scales = 'free') +
   ylab("gas conc (ug/L), flux (ug/L/d)") + 
   xlab('Predictors')+
   ggtitle('Biogeochemical Drivers') +
   theme(legend.position = 'none')
phys <-ggplot(phys_preds, aes(value, gas_value, col = rate)) +
   geom_point() +
   geom_smooth(method = lm)  +
   facet_grid(gas~predictor, scales = 'free') +
   ylab("gas conc (ug/L), flux (ug/L/d)") + 
   xlab('Predictors')+
   ggtitle('Physical Drivers')

png('figures/lm_drivers_individual_linear_regression.png', width = 10, 
    height = 5, res = 300, units = 'in', family = 'cairo')
  ggarrange(bgc, phys, ncol = 2)   
dev.off()

png('figures/lm_bgc_drivers_individual_linear_regression.png', width = 4, 
    height = 5, res = 300, units = 'in', family = 'cairo')
   
dev.off()

# png('figures/lme_drivers_individual_linear_regression.png', width = 8, 
#     height = 5, res = 300, units = 'in', family = 'cairo')
# dat %>%
#   pivot_longer(cols = ends_with(c('.ugL', 'ugld')), names_to = c('gas', 'rate'),
#                                 values_to = 'gas_value', 
#                                 names_pattern = '([A-Z0-9]+).([a-zA-Z_]+)') %>%
#   mutate(value = 
#            case_when(gas == 'DO' & predictor %in%
#                        unique(conc$predictor[conc$gas == "O2"]) ~ value,
#                      gas == 'CO2' & rate == 'ugL' & predictor %in% 
#                        unique(conc$predictor[conc$gas == 'CO2']) ~ value,
#                      gas == 'CO2' & rate == 'flux_ugld' & predictor %in% 
#                        unique(concF$predictor[concF$gas == 'CO2']) ~ value,
#                      gas == 'CH4' & rate == 'ugL' & predictor %in% 
#                        unique(conc$predictor[conc$gas == 'CH4']) ~ value,
#                      gas == 'CH4' & rate == 'flux_ugld' & predictor %in% 
#                        unique(concF$predictor[concF$gas == 'CH4']) ~ value,
#                      gas == 'N2O' & rate == 'ugL' & predictor %in% 
#                        unique(conc$predictor[conc$gas == 'N2O']) ~ value,
#                      gas == 'N2O' & rate == 'flux_ugld' & predictor %in% 
#                        unique(concF$predictor[concF$gas == 'N2O']) ~ value,
#                      TRUE ~ NA_real_),
#          rate = ifelse(rate == 'ugL', 'conc', 'flux')) %>%
#   
#   # filter(rate == 'ugL') 
#  ggplot(aes(value, gas_value, col = rate)) +
#    geom_point() +
#    geom_smooth(method = lm)  +
#    facet_grid(gas~predictor, scales = 'free') +
#    ylab("gas (ug/L)") + 
#    xlab('Predictors')
# dev.off()
# ylim(ylim) 
   geom_text(data = labs, aes(x, y, label = lab), vjust = 1, hjust = 0) +
   ggtitle(paste0(rlab[2], ",  Marginal R2 = ", round(rlab[1], 3)))
for(gg in unique(conc$gas)){
  tmp <- conc %>% 
    filter(gas ==gg) %>%
    arrange(estimate, p_val)
    np <- nrow(tmp)
    preds <- tmp$predictor[1:np]
    rlab <- rsqs %>% filter(gas == gg, flux == T) %>%
      select(R2m, formula)
    
  if(gg == 'O2'){gg = "DO"}
    tt <- dat %>%
      select(gas = any_of(paste0(gg, '.ugL')), predictor, value) %>%
      filter(predictor %in% preds) %>%
      mutate(predictor = factor(predictor, levels = preds)) 
    ylim = c(min(tt$gas), max(tt$gas) + 0.1 * (max(tt$gas) - min(tt$gas)))
    labs = group_by(tt, predictor) %>% 
      summarize(x = min(value, na.rm = T), 
                y = ylim[2])
    for(i in 1:np){
      labs$lab[i] <- tt %>% 
        filter(predictor == preds[i])%>%
        lm_nums()
    }
    p <- ggplot(tt, aes(value, gas)) +
        geom_point() +
        geom_smooth(method = lm, col = 1) +
        facet_wrap(~predictor, nrow = 1, scales = 'free_x') +
        ylab(paste(gg, "(ug/L)")) +
        ylim(ylim) +
        geom_text(data = labs, aes(x, y, label = lab), vjust = 1, hjust = 0) +
        ggtitle(paste0(rlab[2], ",  Marginal R2 = ", round(rlab[1], 3)))
    print(p)
}
              SM_output$K600_daily_Rhat,
              rep('daily', dim(SM_output)[1]))

# DO psat
ann = data.frame()
for(y in 2017:2019){
  m <- summary(lm(NEP ~ DO.obs/DO.sat, data = nhc_fall1[nhc_fall1$year == y,]))
  r = round(m$adj.r.squared, 2)
  p = round(m$coefficients[2,4],2)
  s = round(m$coefficients[2,1],2)
  # mf = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  mf = paste0("r2 = ", r, ",  p =  ", p)
  m <- summary(lm(NEP ~ DO.obs/DO.sat, data = nhc_fall[nhc_fall$year == y,]))
  r = round(m$adj.r.squared, 2)
  p = round(m$coefficients[2,4],2)
  s = round(m$coefficients[2,1],2)
  # my = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  my = paste0("r2 = ", r, ",  p =  ", p)
  ann = bind_rows(ann, data.frame(year = y, mf = mf, my = my))
}
png("figures/NEP_drivers_DO.png", width = 9, height = 3.65, 
    res = 300, units = 'in')
ggplot(nhc_fall, aes(DO.obs/DO.sat, NEP), col = 1) +
  geom_point(size = 2) +
  facet_wrap(.~year) +
  geom_smooth(method = lm, se = F, col = 1) +
  geom_point(data = nhc_fall1, col = fall_col, size = 2) +
  geom_smooth(data = nhc_fall1, method = lm, se =F, col = fall_col) +
  theme_bw() +
  geom_text(data = ann, x = 0.5, y = 5.4, aes(label = mf),
            col = fall_col, hjust = 0) +
  geom_text(data = ann, x = 0.5, y = 5, aes(label = my),
            col = 1, hjust = 0) +
  ylab("") +
  xlab("DO (% saturation)") +
  scale_x_continuous(breaks=c(0,.25,.5,.75,1),
                     labels=c("0","25","50", "75","100"))
dev.off()