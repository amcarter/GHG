# Analyze NHC CO2 - O2 sensor data
# A Carter
# 2020-01-09

library(tidyverse)
library(lubridate)
library(xts)
library(dygraphs)
library(AquaEnv)
library(viridis)
library(streamMetabolizer)
library(ggpubr)

setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")

airpres <- read_csv("data/siteData/interpolatedQ_allsites.csv") %>%
  select(DateTime_UTC, AirPres_kPa, NHC.Q)
dd <- read_csv("data/metabolism/processed/NHC.csv", guess_max = 10000) %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  full_join(airpres, by = "DateTime_UTC") %>%
  select(DateTime_EST, CO2_ppm, DO.obs, DO.sat, CDOM_mV,
         discharge, temp.water, SpecCond_mScm, depth,
         pH, AirTemp_C, AirPres_kPa, Light_lux, NHC.Q, Battery_V) %>%
  arrange(DateTime_EST)

ggplot(dd, aes(Battery_V, CO2_ppm, color = (DateTime_EST))) + 
  ylim(0, 650) +
  geom_point()

dat <- dd %>% 
  # filter((DateTime_EST > ymd_hms("2018-10-01 13:00:00") & 
  #         DateTime_EST < ymd_hms("2019-03-30 00:00:00"))) %>%
  mutate(CO2_ppm = ifelse(CO2_ppm > 9500, NA, CO2_ppm)) %>%
  filter(!is.na(DateTime_EST))

dat$CO2_ppm[(dat$DateTime_EST > ymd_hms("2019-07-21 00:00:00") &
               dat$DateTime_EST < ymd_hms("2019-09-05 15:00:00"))] <- NA 
# Explore CO2 data. Is it reasonable? 
# Note: ppm = uatm = umol/L
hist(dat$CO2_ppm)
dat %>%
  select(CO2_ppm, CDOM_mV, DO.obs) %>%
  xts(order.by = dat$DateTime_EST) %>% 
  dygraph() %>%
  dyRangeSelector() 

# Calculate CO2 saturation:

# Replace this with the actual atmospheric concentration of CO2
CO2_atm <- 450

# calculate deviations from saturation concentration in umol/L
gas <- dat %>% 
  mutate(CO2_ppm = ifelse(CO2_ppm >50000, NA, CO2_ppm),
         CO2.sat_umol = K0_CO2(0, dat$temp.water) * CO2_atm,
         DO_umol = DO.obs * 1000 / 32,
         DO.sat_umol = DO.sat * 1000 / 32, 
         del_CO2 = CO2_ppm - CO2.sat_umol, 
         del_O2 = DO_umol - DO.sat_umol, 
         date = as.Date(DateTime_EST))
png("figures/gas/CO2O2_sensor_NHC2019_full_temp.png", width = 6, height = 6,
    units = "in", res = 300)

ggplot(gas, aes(del_CO2, del_O2, color = log(discharge))) +
  geom_point() +
  geom_abline(slope = -1, intercept = 0) +
  xlab("excess CO2 (umol/L)") +
  ylab("excess O2 (umol/L)") +
  scale_color_gradientn(colors = plasma(8), name = "temp C")
dev.off()
p2 <- ggplot(gas, aes(DateTime_EST, del_O2, color = log(NHC.Q))) +
  geom_point()+
  ylab("excess O2 (umol/L)") +
  scale_color_gradientn(colors = plasma(8), name = "log(Q)")
p3 <- ggplot(gas, aes(DateTime_EST, del_CO2, color = log(NHC.Q))) +
  geom_point()+
  ylab("excess CO2 (umol/L)") +
  scale_color_gradientn(colors = plasma(8), name = "log(Q)")

png("figures/gas/CO2O2_NHC2019_Q.png", width = 6, height = 6,
    units = "in", res = 300)
ggarrange(p3, p2, common.legend = T, ncol = 1)
dev.off()
# add points for grab samples
ghg <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")
ghg_nhc <- ghg %>% 
  mutate(del_CO2 = (CO2.ugL - CO2.sat) / 44, 
         del_O2 = (DO.obs - DO.sat)*1000/32) %>%
#  filter(site == "NHC") %>%
  select(del_CO2, del_O2, site, discharge, watertemp_C, ER, GPP, CH4.ugL) %>%
  filter(site != "MC751")

png("figures/gas/CO2O2_sensor_NHC2019.png", width = 6, height = 6,
    units = "in", res = 300)
ggplot(gas, aes(del_CO2, del_O2, color = temp.water)) +
  geom_point() +
  geom_abline(slope = -1, intercept = 0) +
  scale_color_gradientn(colors = plasma(6)) +
  xlab("excess CO2 (umol/L)") +
  ylab("excess O2 (umol/L)")
  
# geom_point(data = ghg_nhc, color = "black", size = 2)
dev.off()

png("figures/CO2O2_grab_NEP.png", width = 5, height = 4.5,
    units = "in", res = 300)
ggplot(ghg_nhc, aes(del_CO2, del_O2, color = -ER+GPP)) +
  geom_abline(slope = -1, intercept = 0) +
  geom_point(size = 2) +
 scale_color_gradientn(colors = c(plasma(7)[c(4,6)],'forestgreen'),
                       name = 'NEP', na.value = "grey80")+#plasma(6)[1:5]) +
  xlab("excess CO2 (umol/L)") +
  ylab("excess O2 (umol/L)") +
  # theme(legend.position = "bottom")+
  xlim(0, 190)+
  theme_bw()
 ggplot(ghg_nhc, aes(del_CO2, del_O2, color = factor(site))) +
  geom_point(size = 1.5) +
  geom_abline(slope = -1, intercept = 0) +
  geom_smooth(method = lm, se = F) +
 # scale_color_gradientn(colors = plasma(6)[1:5]) +
  xlab("excess CO2 (umol/L)") +
  ylab("excess O2 (umol/L)") +
  theme(legend.position = "bottom")+
  xlim(0, 200)
 
  ggarrange(Q, s) 
dev.off()

ggplot(gas, aes(CO2_ppm- CO2.sat, DO.obs, color = temp.water)) +
  geom_point()
mPULSE package from github

plot(dat$DateTime_EST, dat$NHC.Q, log = "y")
lines(dat$DateTime_EST, dat$discharge)

# The StreamPULSE package is in development and changes frequently!
# If something doesn't work as expected, first try reinstalling.

#library(devtools)
#install_github('streampulse/StreamPULSE', dependencies=TRUE)
