# Prep data frame with variables for NHC ghg analysis
# 2021 01 03

library(tidyverse)
library(lubridate)
library(zoo)

setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")

# load ghg data
ghg <- read_csv("data/gas/NHC_2019-2020_processed_GHGdata.csv") %>%
  rename(watertemp_C.ysi = watertemp_C, 
         airpres_mmHg.ysi = airpres_mmHg, 
         pH.ysi = pH) 

# Combine duplicates examine leaky measures
# To Do: Better quantify measurement error 
# (GC error, how to get obs error from duplicates)
# Account for error where mass/temperature is estimated
ysi <- read_csv("data/water_chemistry/all_nhc_ysi_data.csv") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         DO_pctsat = ifelse(DO_pctsat <= 2, DO_pctsat, DO_pctsat/100),
         DO.sat.ysi = DO_mgL/DO_pctsat) %>%
  dplyr::select(site, date, spc_uScm.ysi = spc_uscm, DO.obs.ysi = DO_mgL, 
         DO.sat.ysi) 

# These dates are grouped so that they are all on the same days longitudinally
ghg <- ghg %>%
  mutate(datetime = round_date(ymd_hms(paste(date, time), tz = 'EST'), 
         unit = "15 minute")) %>%
  dplyr::select(-Notes, -date, -time) %>%
  group_by(site, datetime) %>%
  summarize(across(everything(), mean, na.rm = T)) %>%
  ungroup() %>%
  mutate(date = as.Date(datetime, tz = "EST")) %>%
  left_join(ysi, by = c('site', 'date')) %>%
  mutate(group = case_when(date == as.Date("2019-12-04") ~ "2019-12-03",
                           date == as.Date("2020-01-30") ~ "2020-01-29",
                           TRUE ~ as.character(date)))

# site data
site_dat <- read_csv("data/site_data/NHCsite_metadata.csv") %>%
  slice(c(1:5,7:8)) %>%
  select(site = sitecode, latitude, longitude, CRS, 
         distance_upstream_m = distance_m,
         width_march_m = width_mar_m,
         ws_area_km2 = ws_area.km2,
         slope_wbx, habitat)
  
# instantaneous data
air <- read_csv("data/site_data/NOAA_airpres.csv") %>%
  mutate(datetime = with_tz(DateTime_UTC, tz = 'EST')) %>%
  select(datetime, AirPres_kPa = air_kPa, AirTemp_C = air_temp) 

sites <- site_dat[1:6, 1, drop = T]
# raw <- data.frame()
# for(site in sites){
#   dat <- read_csv(paste0("../nhc_50yl/data/metabolism/processed/", site, ".csv"), 
#                   guess_max = 10000)
#   dat$site <- site
#   raw = bind_rows(raw, dat)
# }
# 
# write_csv(raw, 'data/site_data/metabolism.csv')
raw <- read_csv('data/site_data/metabolism.csv')

raw_dat <- raw %>%
  mutate(spc_uScm = ifelse(!is.na(SpecCond_uScm),
                           SpecCond_uScm, SpecCond_mScm*1000),
         datetime = with_tz(DateTime_UTC, tz = "EST"),
         spc_uScm = case_when(spc_uScm >250 ~ NA_real_,
                              spc_uScm <0.01 ~ NA_real_,
                              TRUE ~ spc_uScm)) %>%
  select(datetime, site, level_m, depth, site, spc_uScm, DO.obs, DO.sat,
         watertemp_C = temp.water, discharge) %>%
  group_by(site) %>%
  mutate(across(-datetime, ~na.approx(., na.rm = F, x = datetime))) %>%
  ungroup()

dat_inst <- ghg %>% 
  left_join(raw_dat, by = c("datetime", "site")) %>%
  left_join(air, by = "datetime") %>%
  mutate(spc_uScm = ifelse(!is.na(spc_uScm), spc_uScm, spc_uScm.ysi),
         DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs.ysi),
         DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat.ysi),
         watertemp_C = ifelse(!is.na(watertemp_C), watertemp_C, 
                              watertemp_C.ysi),
         airpres_kPa = case_when(date == as.Date('2020-02-27') ~ AirPres_kPa,
                              !is.na(airpres_mmHg.ysi) ~ airpres_mmHg.ysi/7.501, 
                              TRUE ~ AirPres_kPa),
         airpres_kPa = (airpres_kPa + AirPres_kPa)/2) %>%
  rename(pH = pH.ysi, airtemp_C = AirTemp_C) %>%
  select(site, date, datetime, group, ends_with(c('ugL')), 
         CH4.sat, CO2.sat, N2O.sat, pH, level_m, 
         discharge, depth, watertemp_C, DO.obs, DO.sat, spc_uScm, 
         airtemp_C, airpres_kPa) 
colnames(dat_inst)  <- c(colnames(dat_inst)[1:11], 
                         paste0(colnames(dat_inst)[12:20], '.inst'))

dd <- dat_inst %>%
  filter(site!='MC751') 
apply(dd, 2, function(x) sum(is.na(x)))
avs <- dd %>%
  filter(!(site %in% c("NHC", "UNHC"))) %>%
  group_by(group) %>%
  summarize(pH.day = mean(pH, na.rm = TRUE),
            spc.day = mean(spc_uScm.inst, na.rm = TRUE),
            DO.obs.day = mean(DO.obs.inst, na.rm = TRUE),
            DO.sat.day = mean(DO.sat.inst, na.rm = TRUE),) 

dd <- dd %>% 
  left_join(site_dat, by = "site") %>%
  left_join(avs, by = 'group') %>%
  mutate(pH = case_when(is.na(pH) ~ pH.day,
                        TRUE ~ pH),
         spc_uScm.inst = case_when(is.na(spc_uScm.inst) ~ spc.day,
                                   TRUE ~ spc_uScm.inst),
         DO.obs.inst = case_when(is.na(DO.obs.inst) ~ DO.obs.day,
                                   TRUE ~ DO.obs.inst),
         DO.sat.inst = case_when(is.na(DO.sat.inst) ~ DO.sat.day,
                                 TRUE ~ DO.sat.inst)) %>%
  group_by(date) %>%
  mutate(discharge.inst = na.approx(discharge.inst, na.rm = F, 
                                    x = ws_area_km2)) %>%
  ungroup() %>%
  select(-ends_with('day'))

dat_inst <- dat_inst %>%
  filter(site == 'MC751') %>%
  left_join(site_dat, by = 'site') %>%
  bind_rows(dd)

# daily data
raw_daily <- raw_dat %>%
  mutate(date = as.Date(datetime, tz = "EST")) %>%
  left_join(air, by = 'datetime') %>%
  rename(airtemp_C = AirTemp_C, airpres_kPa = AirPres_kPa) %>%
  group_by(date, site) %>%
  summarize(across(-datetime, mean, na.rm = T)) %>%
  ungroup() %>%
  select(-level_m)

qq <- left_join(ghg, raw_daily, by = c('site', 'date')) %>%
  select(site, date, discharge) %>%
  pivot_wider(names_from = site, values_from = discharge) %>%
  mutate(NHCQ = NHC) %>%
  pivot_longer(c(-date, -NHCQ), names_to = "site", values_to = "discharge")
Q <- raw_daily %>%
  select(site, date, discharge) %>%
  pivot_wider(names_from = site, values_from = discharge)

for(ss in c('CBP', 'UNHC', 'PM')){
  m <- summary(lm(log(Q[,ss, drop = T])~log(Q$NHC)))$coefficients
  qq <- qq %>%
    mutate(discharge = case_when(!is.na(discharge) ~ discharge,
                                 site == ss ~ exp(m[1,1]) * NHCQ ^ m[2,1]))
}
qq <- select(qq, -NHCQ)
                        
K600 <- readRDS("data/site_data/met_preds_stream_metabolizer.rds")$preds %>%
  filter(era == 'now') %>%
  select(site, date, discharge,starts_with('K600')) %>% 
  distinct()
qq <- bind_rows(K600, qq) %>%
  group_by(site, date) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(across(starts_with('K'), ~ na.approx(., na.rm = F, x = discharge))) %>%
  ungroup()
  
met <- readRDS("data/site_data/met_preds_stream_metabolizer.rds")$filled %>%
  filter(year > 2000) %>% distinct() %>%
  select(site, date, starts_with(c("GPP", "ER"), ignore.case = FALSE)) %>%
  select(-ends_with('cum')) %>%
  # convert to units of O2 from C
  mutate(across(starts_with(c('GPP', 'ER')), ~. * 32/12)) %>%
  left_join(K600[,-3], by = c('site', 'date')) %>%
  group_by(site) %>% 
  mutate(across(-date, ~zoo::na.approx(., na.rm = F, x =  date))) %>% 
  ungroup()

nuts <- read_csv("data/water_chemistry/water_chemistry_2019-2020_compiled.csv") %>%
  select(-sample_name, -time) %>%
  mutate(site = toupper(site))

SP_wchem <- read_csv('data/water_chemistry/StreampulseWQDec2020.csv') %>%
    filter(site %in% c('NHC', 'UNHC')) %>%
    select(site, date, cl_mgl = Cl, so4_mgl = 'SO4 (mg/L)', br_mgl = Br, 
           no3n_mgl = 'NO3-N', na_mgl = 'Na (mg/L)', k_mgl = 'K (mg/L)', 
           mg_mgl = 'Mg (mg/L)', ca_mgl = 'Ca (mg/L)', doc_mgl = 'DOC (mg/L)', 
           tdn_mgl = 'TDN (mg/L)', nh4n_mgl = 'NH4-N (mg/L)', 
           po4p_mgl = 'PO4-P (mg/L)') %>%
    mutate(date = as.Date(date, format = '%m/%d/%Y'),           
           no3n_mgl = case_when(no3n_mgl == "<0.001" ~ "0.0005", 
                                TRUE ~ no3n_mgl),
           no3n_mgl = as.numeric(no3n_mgl)) %>%
    bind_rows(nuts) %>%
    mutate(br_mgl = case_when(br_mgl == "<0.03" ~ "0.015", 
                              br_mgl == "<0.01" ~ "0.005",
                              TRUE ~ br_mgl),
           nh4n_mgl = ifelse(nh4n_mgl == "<0.01", "0.005", nh4n_mgl),
           po4p_mgl = ifelse(po4p_mgl == "<0.01", "0.005",po4p_mgl),
           br_mgl = as.numeric(br_mgl),
           nh4n_mgl = as.numeric(nh4n_mgl),
           po4p_mgl = as.numeric(po4p_mgl),
           date = case_when(site == 'NHC' & date == as.Date('2020-01-29') ~
                              as.Date('2020-01-30'),
                            date == as.Date('2020-03-10') ~ 
                              as.Date('2020-03-11'),
                            TRUE ~ date))
dat_ghg <- dat_inst %>%
  left_join(SP_wchem, by = c("date", "site")) %>%
  left_join(met, by = c("date", "site")) %>%
  left_join(raw_daily, by = c('site', 'date')) %>%
  arrange(date, distance_upstream_m)

colnames(qq) <- c('site', 'date', paste0(colnames(qq)[3:6], '.2'))

dat_ghg <- dat_ghg %>%
  left_join(qq, by = c('site', 'date')) %>%
  mutate(discharge = ifelse(!is.na(discharge), discharge, discharge.2),
         K600 = ifelse(!is.na(K600), K600, K600.2),
         K600_2.5 = ifelse(!is.na(K600_2.5), K600_2.5, K600_2.5.2),
         K600_97.5 = ifelse(!is.na(K600_97.5), K600_97.5, K600_97.5.2)) %>%
  select(-ends_with('.2'))

dat_ghg<- dat_ghg %>% 
  mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs.inst),
         DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat.inst),
         spc_uScm = ifelse(!is.na(spc_uScm), spc_uScm, spc_uScm.inst),
         watertemp_C = ifelse(!is.na(watertemp_C), watertemp_C, 
                              watertemp_C.inst))

date_rng = range(ghg$date)
dat_ghg_predictors <- SP_wchem %>%
  full_join(met, by = c("date", "site")) %>%
  full_join(raw_daily, by = c('site', 'date')) %>%
  left_join(site_dat, by = 'site') %>%
  filter(date <= date_rng[2], date >= date_rng[1]) %>%
  arrange(date, distance_upstream_m)

# filter(dat_ghg, site != 'MC751') %>%
# ggplot(aes(date, ER, col = site)) +
#   geom_point() +
#   geom_line()
# apply( 2, function(x) sum(is.na(x)))
  
write_csv(dat_ghg, 'data/ghg_complete_drivers_dataframe.csv')
write_csv(dat_ghg_predictors, 'data/ghg_filled_drivers_dataframe.csv')
  
# calculate gas flux ####
library(marelac)
# K_gas = K600(Sc/600)^-0.5
K <-  data.frame(K600 = dat_ghg$K600,
                 gas_schmidt(dat_ghg$watertemp_C.inst, 
                             c("CH4","CO2", "N2O", "O2"))) %>%
  mutate(across( -K600, ~(./600)^(-0.5)*K600, 
                .names = "K_{col}")) %>%
  select(starts_with('K_')) %>%
  as_tibble() %>%
  bind_cols(dat_ghg) %>%
  mutate(across(ends_with('ugL'), ~ifelse(.<0, 0, .)))%>%
  mutate(CH4.flux_ugld = (CH4.ugL - CH4.sat) * K_CH4,
         CO2.flux_ugld = (CO2.ugL - CO2.sat) * K_CO2,
         N2O.flux_ugld = (N2O.ugL - N2O.sat) * K_N2O,
         O2.flux_ugld = (DO.obs - DO.sat) * K_O2 * 1000) 

write_csv(K, "data/ghg_flux_complete_drivers_dataframe.csv")
dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")

# pare down to have fewer NA's
d2 <- dat %>%
  filter(site != "MC751", 
         !is.na(K600),
         !is.na(CH4.ugL),
         !is.na(GPP)) %>%
  select(-pH, -level_m.inst, -na_mgl, -mg_mgl, -ca_mgl, -br_mgl,
         -starts_with("K_"), -CRS, -latitude, -longitude, -po4p_mgl, -cl_mgl, 
         -datetime, -so4_mgl, -spc_uScm, -starts_with('air')) #%>%
  # filter(!is.na(GPP)) 


# ggplot(d2, aes(discharge, discharge.inst, col = site)) +
#   geom_line()+
#   geom_point()
#   facet_wrap(.~site)
apply(d2, 2, function(x) sum(is.na(x)))
write_csv(d2, "data/ghg_flux_complete_drivers_dataframe_noNAs.csv")
