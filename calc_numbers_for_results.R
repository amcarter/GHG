#summarize numbers for results
library(tidyverse)
library(lubridate)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")

dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")%>%
  filter(!is.na(datetime),
         site != "MC751")
dvs <- read_csv("data/ghg_filled_drivers_dataframe.csv")


# magnitudes ####
dat$site <- factor(dat$site, levels = c("NHC", "PM", "CBP", "WB", "WBP","UNHC"))
dat <- dat %>% mutate(across(ends_with('flux_ugld'), ~.*depth)) %>%
  rename_with(ends_with("flux_ugld"),.fn = ~gsub("_ugld", "_mgm2d", .) )

dat %>% dplyr::select(ends_with(c('.obs','ugL', 'mgm2d'))) %>%
  summary()

# anova ####
anova(lm(CO2.ugL ~ site, data = dat ))
anova(lm(CH4.ugL ~ site, data = dat ))
anova(lm(N2O.ugL ~ site, data = dat ))
anova(lm(DO.obs ~ site, data = dat ))

anova(lm(CO2.flux_mgm2d ~ site, data = dat ))
anova(lm(CH4.flux_mgm2d ~ site, data = dat ))
anova(lm(N2O.flux_mgm2d ~ site, data = dat ))
anova(lm(O2.flux_mgm2d ~ site, data = dat ))
# seasonal patterns
unique(dat$group)
dat %>% dplyr::select(c(group, DO.obs)) %>%
  filter(group %in% as.Date(c('2019-11-11','2020-01-29', '2020-03-20'))) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, CO2.ugL)) %>%
  filter(group %in% as.Date(c('2019-11-11','2019-12-03',
                              '2019-12-12','2020-03-11'))) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, CO2.ugL)) %>%
  filter(group %in% as.Date(c('2019-12-03',
                              '2019-12-12'))) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))


dat %>% dplyr::select(c(group, CH4.ugL)) %>%
  filter(group %in% as.Date(c('2019-11-11','2019-12-03',
                              '2019-12-12','2020-03-20'))) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, CH4.ugL)) %>%
  filter(group %in% as.Date(c('2019-12-03',
                              '2019-12-12','2020-01-05', '2020-01-29'))) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, N2O.ugL)) %>%
  # filter(group %in% as.Date(c('2019-11-11','2019-12-03',
  #                             '2019-12-12','2020-03-20'))) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
# fluxes ####
dat %>% dplyr::select(c(group, O2.flux_mgm2d)) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, O2.flux_mgm2d)) %>%
  filter(group %in% as.Date(c('2019-11-11',
                              '2019-11-20'))) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))

dat %>% dplyr::select(c(group, CO2.flux_mgm2d)) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, CO2.flux_mgm2d)) %>%
  filter(group %in% as.Date(c('2019-11-11',
                              '2019-11-20',
                              '2019-11-26'))) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))

dat %>% dplyr::select(c(group, CH4.flux_mgm2d)) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))
dat %>% dplyr::select(c(group, CH4.flux_mgm2d)) %>%
  filter(group %in% as.Date(c('2019-11-11',
                              '2019-11-20'))) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))

dat %>% dplyr::select(c(group, N2O.flux_mgm2d)) %>%
  group_by(group) %>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                            sd = ~sd(., na.rm = T),
                            n = ~n()))

ins <- read_csv('data/fraction_of_instream_production_CO2_and_CH4CO2ratios.csv')
ss <- ins %>%
  select(site, date, CO2_flux, NEP, Extra, instr) %>%
  group_by(site) %>%
  select( -date)%>%
  summarize_all(.funs = list(mean = ~mean(., na.rm = T),
                sd = ~sd(., na.rm = T)))

summary(aov(instr ~CO2_flux, data = ins))

filter(CO2_flux>0) %>%
  mutate(instr = ifelse(instr>1, 1, instr)) %>%
  summarize_all(sum, na.rm = T)
  summary()
<- dat %>%
  filter(!is.na(datetime),
         site != "MC751") %>%
  # mutate(date = as.Date(group))%>%
  select(datetime, date, site,  habitat, distance_upstream_m, watertemp_C,
         depth, GPP, ER, discharge, DO.obs, no3n_mgl, ends_with('ugld')) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas",
               names_pattern ='([0-9A-Z]+).', values_to = "flux_ugld") %>%
  left_join(gas, by = c("date",'datetime', "site", "gas")) %>%
  mutate(gas = factor(gas, levels = c('CO2','O2',  'CH4', 'N2O')),
         flux_mgm2d = flux_ugld * depth)%>%
  arrange(date, distance_upstream_m)

