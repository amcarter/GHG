# explore driver variables
library(tidyverse)
library(lubridate)
library(car)

setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")
# setwd('C://Users/adelv/Dropbox/Duke/NHC/Alice')

dvs <- read_csv("data/ghg_filled_drivers_dataframe.csv") %>%
  filter(site !='MC751') %>%
  mutate(site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC'))) 

dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv") %>%
  filter(site !='MC751') %>%
  mutate(date = as.Date(group),
         site = factor(site, levels = c('UNHC', 'WBP','WB','CBP','PM','NHC'))) 


summary(lm(DO.obs~site-1, data = dat))
TukeyHSD(aov(K600~site, data = dat))
pairwise.t.test(dat$K600, dat$site)
unique(dat$date)
dat %>% filter(!(site %in% c('UNHC'))) %>% 
  summarize(mean = mean(no3n_mgl, na.rm = T), 
            sd = sd(no3n_mgl, na.rm = T))
dat %>% filter(date > as.Date('2020-03-01')) %>% 
  summarize(mean = mean(no3n_mgl, na.rm = T), 
            sd = sd(no3n_mgl, na.rm = T),
            n = n())
# d2 <- read_csv("data/ghg_flux_complete_drivers_dataframe_noNAs.csv")
  
ggplot(dat, aes(site, DO.obs, fill = site)) +
  geom_boxplot()

ggplot(dat, aes(factor(date, DO.obs))+#, col = site)) +
  geom_boxplot()
ggplot(dat, aes(date, DO.obs, col = site)) +
  geom_line()

# summarize discrete drivers ####
#DO ####
dat %>% filter((site %in% c('NHC', 'PM'))) %>% 
  summarize(mean = mean(DO.obs, na.rm = T), 
            sd = sd(DO.obs, na.rm = T), n = n())
dat %>% filter(date > as.Date('2020-03-01')) %>% 
  summarize(mean = mean(no3n_mgl, na.rm = T), 
            sd = sd(no3n_mgl, na.rm = T),
            n = n())

# summarize continuous drivers ####
# depth ####
dvs %>% filter((site %in% c('WB', 'WBP'))) %>% 
  summarize(mean = mean(depth, na.rm = T), 
            sd = sd(depth, na.rm = T))
# discharge ####
dvs %>% filter(date >  as.Date('2020-03-10')) %>% 
  ggplot(aes(date, watertemp_C, col = site)) +
  geom_line()# + ylim(0,3)
  # group_by(site) %>%
  summarize(mean = mean(discharge, na.rm = T), 
            sd = sd(discharge, na.rm = T),
            n = n())
# bottom of reach in mid Nov
dvs %>% filter(site == 'NHC',
                 date < as.Date('2019-11-24')) %>% 
  summarize(mean = mean(discharge, na.rm = T), 
            sd = sd(discharge, na.rm = T),
            n = n())
# late Nov - mid Feb with 3 storms
dvs %>% filter(site == 'NHC',
               date > as.Date('2019-11-23') &
                 date < as.Date('2020-02-16'),
               discharge < 3) %>% 
  summarize(mean = mean(discharge, na.rm = T), 
            sd = sd(discharge, na.rm = T),
            n = n())
# late Jan - mid Feb two large and one small storms
dvs %>% filter(site == 'NHC',
               date %in% as.Date(c('2020-01-04', '2020-01-25', 
                                   '2020-02-06'))) %>% 
  dplyr::select(date, discharge)
# late Nov - mid Feb with 3 storms
dvs %>% filter(site == 'NHC',
               date >= as.Date('2020-02-16')) %>% 
  summarize(mean = mean(discharge, na.rm = T), 
            sd = sd(discharge, na.rm = T),
            n = n())

# temperature ####

# Nov Dec
dvs %>% filter(date < as.Date('2019-12-27')) %>% 
  summarize(mean = mean(watertemp_C, na.rm = T), 
            sd = sd(watertemp_C, na.rm = T),
            n = n())
# Jan range
dvs %>% filter(date > as.Date('2019-12-26')&
                 date < as.Date('2020-02-01')) %>% 
  dplyr::select(watertemp_C) %>% summary()
# feb mar
dvs %>% filter(date > as.Date('2020-01-31') &
                 date < as.Date('2020-03-10'),
               watertemp_C >5) %>% 
  summarize(mean = mean(watertemp_C, na.rm = T), 
            sd = sd(watertemp_C, na.rm = T),
            n = n())
#latemar
dvs %>% filter(date > as.Date('2020-03-10')) %>% 
  summarize(mean = mean(watertemp_C, na.rm = T), 
            sd = sd(watertemp_C, na.rm = T),
            n = n())
