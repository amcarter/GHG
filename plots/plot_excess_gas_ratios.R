
library(tidyverse)
library(lubridate)
library(xts)
library(dygraphs)
library(AquaEnv)
library(viridis)
library(streamMetabolizer)
library(ggpubr)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")
ghg <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
ghg_nhc <- ghg %>%
  mutate(del_CO2 = (CO2.ugL - CO2.sat) / 44,
         del_O2 = (DO.obs - DO.sat)*1000/32,
         del_CH4 = (CH4.ugL - CH4.sat)/16,
         del_N2O = (N2O.ugL - N2O.sat)/16
         ) %>%
  #  filter(site == "NHC") %>%
  # select(del_CO2, del_O2, del_CH4, site, discharge, watertemp_C, ER, GPP,
  #        DO.obs, doc_mgl, CH4.ugL) %>%
  filter(site != "MC751")
slope <- read_csv('data/sites_slope_comparison.csv') %>%
  select(site, slope_mm)

ghg_nhc <- ghg_nhc %>% left_join(slope) %>%
  mutate(wrt_days = depth *width_march_m*1000/discharge/60/60/24,
         CH4_rtd = (CH4.ugL * width_march_m * depth * 1000)/
           (CH4.ugL *discharge * 60*60 *24 +
                       CH4.flux_ugld * depth * width_march_m * 1000),
         CO2_rtd = (CO2.ugL * width_march_m * depth * 1000)/
           (CO2.ugL *discharge * 60*60 *24 +
                       CO2.flux_ugld * depth * width_march_m * 1000),
         N2O_rtd = (N2O.ugL * width_march_m * depth * 1000)/
           (N2O.ugL *discharge * 60*60 *24 +
                       N2O.flux_ugld * depth * width_march_m * 1000),
         velocity = discharge/depth/width_march_m,
         CH4_turnover_m = 3 * velocity/ K_CH4 * 60* 60 * 24,
         CO2_turnover_m = 3 * velocity/ K_CO2 * 60* 60 * 24,
         N2O_turnover_m = 3 * velocity/ K_N2O * 60* 60 * 24,
         O2_turnover_m = 3 * velocity/ K_O2 * 60* 60 * 24,
         NEP = GPP - ER
         )
png('figures/CO2_vs_WRT_plot.png', width = 5, height = 4, units = 'in',
    family = 'cairo', res = 300)
ghg_nhc %>%
  mutate(site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC')),
         slope_F = ifelse(slope_mm >0.003, "H", "L"))%>%
  mutate(CO2_turnover_d = CO2_turnover_m * width_march_m * depth/discharge/60/60/24) %>%
  # pivot_longer(ends_with('rtd'), names_to = 'gas', values_to = 'restime',
  #              names_pattern = '([A-Z,0-9]+)_rtd')%>%
ggplot(aes((wrt_days),  CO2_turnover_d, color = site)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2)+
  xlab("Water Residence Time (days)") +
  ylab("CO2 Residence Time (days)") +
  ylim(0, 4)+
  theme_bw()
dev.off()
ghg_nhc %>% pivot_longer(cols = ends_with('turnover_m'),
                         values_to = "turnover_m", names_to = 'gas')%>%
  mutate(turnover_days = turnover_m * width_march_m * depth/discharge/60/60/24) %>%
ggplot(aes((wrt_days), (turnover_days), col  = site)) +
  geom_point() +   geom_smooth(method = lm, se = F)+
  geom_abline(slope = 1, intercept = 0)

# Figure 5 excess CO2CH4 ####
png('figures/excess_O2_CH4_CO2_plot.png', width = 8, height = 4, units = 'in',
    family = 'cairo', res = 300)
  CH <- ggplot(ghg_nhc, aes(del_CO2, del_CH4, color =log10(discharge))) +
    geom_smooth(method = lm, se = F, col = '#999999')+
    geom_point(size = 2) +
    scale_color_gradientn(colors = plasma(8)[1:7],
                          name = 'Discharge (m3/s)',
                          na.value = "grey80",
                          breaks = c(-1,-.30103, .30103),
                          labels = c(0.1,0.5,2)) +
    xlab("CO2 departure (umol/L)") +
    ylab("CH4 departure (umol/L)") +
    theme_bw()+
    theme(legend.key.size = unit(0.5,'cm'),
          legend.title = element_text(size = 8.5),
          legend.position = c(.4, .9),
          legend.direction = 'horizontal')

  DO <- ggplot(ghg_nhc, aes(del_CO2, del_O2, color = GPP-ER)) +
    geom_smooth(method = lm, se = F, col = '#999999')+
    geom_point(size = 2) +
    geom_abline(slope = -1.25, intercept = 0, lty = 2, lwd = .5) +
    scale_color_gradientn(colors = c(plasma(7)[c(4:6)],'forestgreen'),
                          name = 'NEP (gO2/m2/d)', na.value = "grey80")+#plasma(6)[1:5]) +
    xlab("CO2 departure (umol/L)") +
    ylab("O2 departure (umol/L)") +
    theme_bw()+
    theme(legend.key.size = unit(0.5,'cm'),
          legend.title = element_text(size = 8.5),
      legend.position = c(.6, .9),
      legend.direction = 'horizontal')

  ggarrange( DO,CH, ncol = 2)

dev.off()

summary(lm(del_O2 ~ del_CO2, data = ghg_nhc))
ghg_nhc$delll <- ghg_nhc$del_CH4/ghg_nhc$del_CO2
ghg_nhc%>% group_by(group) %>%
    summarize(m = mean(delll))
sd(ghg_nhc$delll)
png("figures/CO2O2_grab_NEP.png", width = 5, height = 3.5,
    units = "in", res = 300)
ggplot(ghg_nhc, aes(del_CO2, del_O2, color = factor(site))) +
  geom_point(size = 1.5) +
  geom_abline(slope = -1, intercept = 0) +
  geom_smooth(method = lm, se = F) +
  # scale_color_gradientn(colors = plasma(6)[1:5]) +
  xlab("excess CO2 (umol/L)") +
  ylab("excess O2 (umol/L)") +
  theme(legend.position = "bottom")+
  xlim(0, 200)

# ggarrange(Q, s)
dev.off()


# expected CO2 based on aerobic respiration

exp <- ghg_nhc %>%
  # select(starts_with('del'),ER, GPP, site, date)%>%
  mutate(exp_CO2 = -del_O2 /1.25,
         extra_CO2 = del_CO2 - exp_CO2,
         site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC')))
exp2 <- ghg_nhc %>%
  # select(starts_with('del'),ER, GPP, site, date)%>%
  mutate(exp_CO2 = -DO.obs*1000/32 /1.25,
         extra_CO2 = CO2.ugL/44 - exp_CO2,
         site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC')))
ggplot(exp2, aes(CH4.ugL, CH4.ugL/16/(extra_CO2), col = date)) +
  geom_point(size = 2) #+ geom_smooth(se = F)
  geom_point(aes(y = del_CH4/extra_CO2), col = 'black') + geom_smooth(col = 'black')

summary(lm(del_CH4 ~ extra_CO2, data = exp))
summary(aov(K600 ~ site, data = exp))
summary(exp$del_CH4/exp$extra_CO2)

png('figures/excess_CO2Met_CH4_plot.png', width = 8, height = 4, units = 'in',
    family = 'cairo', res = 300)

DO <- ggplot(exp, aes(extra_CO2, del_CH4, col = ER)) +
  geom_point(size = 2)+
  scale_color_gradientn(colors = c(plasma(4)), name = 'DO (mg/L)')+
  geom_smooth(method = lm, se = F, col = '#999999') +
  xlab("CO2 departure in excess of NEP (umol/L)") +
  ylab("CH4 departure (umol/L)") +
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'),
        legend.title = element_text(size = 8.5),
        legend.position = c(.27, .92),
        legend.direction = 'horizontal')
ss <- ggplot(exp, aes(extra_CO2, del_CH4, col = site)) +
  geom_point(size = 2)+
  # scale_color_gradientn(colors = c(plasma(5)))+
  geom_smooth(method = lm, se = F) +
  xlab("CO2 departure in excess of NEP (umol/L)") +
  ylab("CH4 departure (umol/L)") +
  theme_bw()+
  theme(legend.key.size = unit(0.4,'cm'),
        legend.title = element_text(size = 8.5),
        legend.position = c(.13, .8),
        legend.direction = 'vertical')

ggarrange(DO, ss,  ncol = 2)

dev.off()


ggplot(exp, aes(DO.obs, del_CH4/del_CO2)) +
  geom_point(col =plasma(4)[3]) +
  geom_smooth(aes(y = del_CH4/del_CO2),method = lm, col =plasma(4)[3], se = F)+
  geom_point(aes(y = del_CH4/extra_CO2), col =plasma(4)[1])+
  geom_smooth(aes(y = del_CH4/extra_CO2),method = lm,col =plasma(4)[1], se = F) +
  ylab('CH4:CO2 (orange),  CH4:extra CO2 (purple)') +
  xlab("DO (mg/L)")+
  theme_bw()

summary(lm(del_O2~NEP, data = exp))
points(exp$del_CO2, exp$del_O2, col = 1, pch = 19)


# CO2 from instream production ####
tmp <- ghg_nhc %>%
  mutate(CO2_flux = (CO2.flux_ugld * depth)/1000,
         NEP_CO2 = ((ER - GPP)*44/32/1.25),
         NEP = ifelse(NEP_CO2>0, NEP_CO2, 0),
         NEP_CO2_high = ((ER - GPP)*44/32 *1),
         NEP_high = ifelse(NEP_CO2_high>0, NEP_CO2_high, 0),
         NEP_CO2_low = ((ER - GPP)*44/32 *0.6),
         NEP_low = ifelse(NEP_CO2_low>0, NEP_CO2_low, 0),
         Extra = case_when(CO2_flux < 0~ NA_real_,
                           CO2_flux < NEP ~ 0,
                           TRUE ~ CO2_flux- NEP),
         Extra_high = case_when(CO2_flux < 0~ NA_real_,
                           CO2_flux < NEP_high ~ 0,
                           TRUE ~ CO2_flux- NEP_high),
         Extra_low = case_when(CO2_flux < 0~ NA_real_,
                           CO2_flux < NEP_low ~ 0,
                           TRUE ~ CO2_flux- NEP_low),
         instr = NEP/CO2_flux,
         # instr = case_when(instr < 0 ~ 0,
         #                   instr > 1 ~ 1,
         #                   TRUE ~instr),
         instr_high = NEP_high/CO2_flux,
         # instr_high = case_when(instr_high < 0 ~ 0,
         #                   instr_high > 1 ~ 1,
         #                   TRUE ~instr_high),
         instr_low = NEP_low/CO2_flux,
         # instr_low = case_when(instr_low < 0 ~ 0,
         #                   instr_low > 1 ~ 1,
         #                   TRUE ~instr_low),
         date = as.Date(group),
         CO2.umol = CO2.ugL/44,
         CH4CO2 = (CH4.ugL/16)/(CO2.umol),
         CO2.umol_NEP = -del_O2 / 1.25,
         CO2.umol_extra = ifelse(CO2.umol_NEP > CO2.umol, 0, CO2.umol - CO2.umol_NEP),
         CH4CO2_extra = (CH4.ugL/16)/(CO2.umol_extra),
         site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC')))
  # select(site, date, instr, CO2_flux, discharge, NEP_CO2, Extra, NEP, CH4CO2)
write_csv(tmp, 'data/fraction_of_instream_production_CO2_and_CH4CO2ratios.csv')

summary(lm(del_CO2 ~ NEP_CO2, data = tmp))
summary(lm(del_O2 ~ NEP_CO2, data = tmp))
11.039/14.58

plot(tmp$date, tmp$CH4CO2, pch = 19)
points(tmp$date, tmp$CH4CO2_extra)

tmp %>%
  pivot_longer(cols = starts_with("CH4CO2"), names_to = "CH4CO2", values_to = "CH4CO2ratio") %>%
  ggplot(aes(date, CH4CO2ratio, col = CH4CO2)) +
  geom_point(size = 2) + geom_smooth(se = F) +
  theme_bw()
tmp %>%
  pivot_longer(cols = any_of(c("NEP", "Extra")), values_to = "gm2d",
               names_to = "category") %>%
  ggplot(aes(date, gm2d, col = category)) +
  geom_point(size = 2) + geom_smooth(se = F) +
  theme_bw()

png("figures/CO2flux_from_instream_production.png", width = 4.5, height = 3.5,
    res = 300, family = 'cairo', units = 'in')

dd <- tmp %>%
  pivot_longer(cols = any_of(c("NEP", "Extra")), values_to = "gm2d",
               names_to = "category") %>%
  ggplot(aes(date, gm2d, fill = category) )+
  geom_bar(stat = 'identity') + #ylim(-.1,1) +
  ylab('CO2 flux (g/m2/d)') +
  xlab('')+
  guides(fill = F)+
  scale_fill_manual(values = c('grey70', alpha('forestgreen', .8)))+
  theme_bw()+
  theme(plot.margin = unit(c(.1,.2,.2,.2), "cm"))
ss <- tmp %>%
  pivot_longer(cols = any_of(c("NEP", "Extra")), values_to = "gm2d",
               names_to = "category") %>%
  ggplot(aes(site, gm2d, fill = category) )+
  geom_bar(stat = 'identity', width = .5) + #ylim(-.1,1) +
  ylab('CO2 flux (g/m2/d)') +
  scale_fill_manual(values = c('grey70', alpha('forestgreen', .8)))+
  theme_bw()+
  theme(plot.margin = unit(c(.2,3,.1,.2), "cm"))
ggarrange(ss, dd, ncol = 1)

dev.off()

tmp %>% group_by(site) %>%
  summarize(across(starts_with('instr'),  mean, na.rm = T))
tmp %>% group_by(date) %>%
  summarize(across(starts_with('instr'),  mean, na.rm = T))

tmp %>% select(site, date, starts_with(c('NEP_CO2', 'Extra'))) %>%
    # filter(site == "UNHC") %>%
    arrange(Extra) %>% filter(NEP_CO2 <0)

tmp %>% filter(NEP_CO2_low > 0) %>% select(NEP_CO2_low) %>% summarize(s = sum(NEP_CO2_low))
tmp %>% filter(NEP_CO2_high > 0) %>% select(NEP_CO2_high) %>% summarize(s = sum(NEP_CO2_high))
summary(tmp)
plot(tmp$NEP_CO2, tmp$CH4CO2)
ggplot(tmp, aes(discharge, instr, col = date)) + geom_point(size = 3)

summary(lm(instr~site, data = tmp))

summary(tmp)

