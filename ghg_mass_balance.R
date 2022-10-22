# Initial GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(ggpubr)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")

dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv") %>%
  mutate(distance_m = 8450 - distance_upstream_m,
         datetime = with_tz(datetime, tz = "EST")) %>%
  filter(site != 'MC751')
d2 <- read_csv("data/ghg_flux_complete_drivers_dataframe_noNAs.csv")
Q <- read_csv('data/ghg_filled_drivers_dataframe.csv') %>%
    filter(site == "NHC") %>%
    # filter(site %in%c("NHC", 'UNHC')) %>%
    mutate(discharge_m3d = discharge * 24*60*60) %>%
    select(date, discharge_m3d)
Q2 <- read_csv('data/ghg_filled_drivers_dataframe.csv') %>%
    filter(site %in%c("NHC", 'UNHC')) %>%
    mutate(discharge_m3d = discharge * 24*60*60) %>%
    select(date, site, discharge_m3d) %>%
    pivot_wider(id_cols = date, names_from = site, values_from = discharge_m3d)
ggplot(dat, aes(ws_area_km2, distance_m)) +
  geom_point()
# Mass Balance Calculations ####
# F_out = F_in + F_met_prod - F_evasion + F_missing

#filter dataframe
gas <- dat %>%
  mutate(O2.ugL = DO.obs.inst * 1000) %>%
  select(date, site, CH4.ugL, CO2.ugL, N2O.ugL, O2.ugL) %>%
  pivot_longer(ends_with("ugL"), names_to = "gas",
               names_pattern = '([A-Z0-9]+).ugL',
               values_to = "concentration_ugl")
sat <- dat %>%
  mutate(O2.sat = DO.sat * 1000) %>%
  select(date, site, ends_with("sat"), -DO.sat) %>%
  pivot_longer(ends_with("sat"), names_to = "gas",
               names_pattern = '([A-Z0-9]+).sat',
               values_to = "sat_ugl")

k <- dat %>%
  select(date, site, starts_with("K_", ignore.case = F)) %>%
  pivot_longer(starts_with("K_"), names_to = "gas",
               names_prefix = 'K_',
               values_to = "K_gas_d")

mb <- dat %>%
  select(date, group, site, distance_m, ws_area_km2,
         watertemp_C, airpres_kPa, airtemp_C, GPP, ER,
         ends_with("ugld")) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas", values_to = "flux_ugld",
               names_pattern = '([A-Z0-9]+).') %>%
  left_join(gas, by = c("date", "site", "gas")) %>%
  left_join(k, by = c("date", "site", "gas")) %>%
  left_join(sat, by = c("date", "site", "gas")) %>%
  mutate(site = factor(site, levels = c("UNHC", "WBP", "WB", "CBP", "PM","NHC")),
         reach_no = as.numeric(site)-1)%>%
  arrange(date, distance_m)

sites <- dat %>%
  select(site, distance_m, width_march_m) %>%
  group_by(site) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  arrange(distance_m)
sites <- sites %>%
  mutate(length_m = c(NA, diff(sites$distance_m)),
         sa_m2 = length_m * width_march_m) %>%
  select(site, length_m, sa_m2) %>%
  full_join(dat) %>%
  select(site, date, distance_m, discharge, depth_m = depth,
         width_m = width_march_m, length_m, sa_m2)# %>%
  # mutate(discharge = ifelse(site == 'NHC' &
  #                             date %in% as.Date(c("2019-11-20",
  #                                                 "2020-01-29", "2020-02-12")),
  #                           NA,discharge))

mb <- left_join(mb, sites, by = c('site', 'date', 'distance_m')) %>%
  mutate(date = as.Date(group)) %>%
  select(-group)
# F_in is the mass flux across the input point:
# F_in = pGas (mg/m3) * discharge (m3/s)
# first, CO2

#
# dd <- mb %>%
#   filter(date == dates[1]) %>%
#   # arrange(distance_m) %>%
#   mutate(discharge_m3s = na.approx(discharge_m3s, ws_area_km2, na.rm = T),
#          depth = na.approx(depth, distance_m, na.rm = F),
#          GPP = na.approx(GPP, distance_m, na.rm = F),
#          ER = na.approx(ER, distance_m, na.rm = F),
#          CH4.ugL = na.approx(CH4.ugL, distance_m, na.rm = F),
#          CO2.ugL = na.approx(CO2.ugL, distance_m, na.rm = F),
#          N2O.ugL = na.approx(N2O.ugL, distance_m, na.rm = F),
#          CH4.flux_ugld = na.approx(CH4.flux_ugld, distance_m, na.rm = F),
#          CO2.flux_ugld = na.approx(CO2.flux_ugld, distance_m, na.rm = F),
#          N2O.flux_ugld = na.approx(N2O.flux_ugld, distance_m, na.rm = F))
#

gO2_to_gCO2 = function(x, PQ){
  x = x * (1 / (2 * 15.9994)) * (1 / PQ) * (12.0107 + (2 * 15.9994))
  return(x)
}

dates <- unique(mb$date)#[c(-7,-8)]
flux <- data.frame()

for(i in 1:length(dates)){
  dd_out <- mb %>%
    filter(date == dates[i], reach_no > 0) %>%
    mutate(Fout_gd = discharge * concentration_ugl *60*60*24/1000,
           Met_gO2m2d = (GPP-ER),
           K_d = K_gas_d) %>%
    select(date, site, reach_no, gas, discharge, distance_m, depth_m, length_m,
           ws_area_km2, sa_m2, c_out_ugl = concentration_ugl,
           csat_out_ugl = sat_ugl,Fout_gd,
           Met_gO2m2d, K_d)
  dd_in <- mb %>%
    filter(date == dates[i], reach_no < 5) %>%
    mutate(reach_no = reach_no + 1,
           Fin_gd = discharge * concentration_ugl *60*60*24/1000) %>%
    select(reach_no, gas, Fin_gd, c_in_ugl = concentration_ugl,
           csat_in_ugl = sat_ugl, discharge_in = discharge)


  dd <- left_join(dd_out, dd_in, by = c('reach_no', 'gas')) %>%
    mutate(F_met_gO2d = Met_gO2m2d * sa_m2,
           F_met_gd = case_when(gas =='O2'~ F_met_gO2d,
                                gas == 'CO2' ~
                                  -gO2_to_gCO2(F_met_gO2d, 1.25),
                                TRUE ~ 0),
           F_atm_gd = -K_d * ((c_in_ugl + c_out_ugl)/2 -
                               (csat_in_ugl + csat_out_ugl)/2) *
                      sa_m2 * depth_m/1000,
           # F_atm_gd = ifelse(site == 'WB', F_atm_gd *2, F_atm_gd),
           F_delta_gd = Fout_gd - Fin_gd,
           F_missing_gd = F_delta_gd - F_met_gd - F_atm_gd,
           delta_Q_m3d = (discharge - discharge_in) *60*60*24) %>%
    select(-F_met_gO2d) %>%
    mutate(across(starts_with('F_'), ~./sa_m2, .names = '{col}2' ),
           delta_Q_md = delta_Q_m3d/sa_m2) %>%
    rename_with(.cols = ends_with('_gd2'),
                .fn = function(x) sub('gd2', 'gm2d', x = x)) %>%
    select(date, site, distance_m, ws_area_km2, reach_no, length_m, discharge, sa_m2,
           gas, starts_with(c('c', 'F', 'delta')))


  flux <- bind_rows(flux, dd)
}

flux <- flux %>%
  mutate(gas = factor(gas, levels = c('CO2','O2', 'CH4', 'N2O')))
f_in <- flux %>%
    filter(site == 'WBP') %>%
    select(date, gas, Fin_gd)

f_inout <- flux %>%
    filter(site == 'NHC') %>%
    select(date, gas, Fout_gd) %>%
    left_join(f_in, by = c('date', 'gas'))

ff <- flux %>% group_by(date, gas) %>%
    select(date, ends_with(c('_gd', '_m3d')), sa_m2, length_m) %>%
    select(-Fin_gd, -Fout_gd) %>%
    summarize_all(.funs = function(x) sum(x, na.rm = T)) %>%
    left_join(Q, by = 'date') %>%
    filter(gas != 'O2') %>%
    left_join(f_inout, by = c('date', 'gas'))

subcatchement_area <- dat$ws_area_km2[dat$site == 'NHC'][1] -
    dat$ws_area_km2[dat$site == 'UNHC'][1]

ff$gw_flux_m3m2d =  ff$delta_Q_m3d/ff$sa_m2

f <- select(ff, date, gw_flux_m3m2d) %>%
    group_by(date) %>%
    summarize(gw_flux_m3m2d = mean(gw_flux_m3m2d))

write_csv(f,'data/gw_fluxes.csv')
Q2 <- Q2 %>% left_join(select(ff, date, sa_m2)) %>%
    mutate(sa_m2 = zoo::na.approx(sa_m2, x = date),
           gw_flux = (NHC - UNHC)/sa_m2)

Q2 %>% group_by(date) %>%
    summarize(gw_flux = mean(gw_flux, na.rm = T)) %>%
    filter(date %in% unique(ff$date))
tiff('figures/final/gw_flux_distribution.tiff',
     height = 4, width = 5, units = 'in', res = 300)
    par(mar = c(3,4.5,1,2),
        mfrow = c(2,1))
    plot(density(ff$gw_flux_m3m2d, na.rm = T),
         xlab = '', ylab = '', cex.axis = 0.8, main = '')
    mtext(expression(paste('Groundwater Flux  (', m^3, m^-2, d^-1, ')')),1,2.1,
          cex = 0.85)
    mtext('Density', 2, 2.1, cex = 0.85)
    points(ff$gw_flux_m3m2d, rep(0.03, nrow(ff)), pch = 17, col = 'brown3')
    legend('topleft', legend = 'sample days', pch = 17, col = 'brown3', bty = 'n')

    # plot(Q2$date, Q2$gw_flux, type = 'l',
    plot(ff$date, ff$gw_flux_m3m2d, type = 'b', pch = 20,
         xaxt = 'n', xlab = 'Date',
         ylab = '', cex.axis = 0.8)
    abline(h = 0, col = 'grey70')
    points(ff$date, ff$gw_flux_m3m2d, pch = 17, col = 'brown3')
    mtext(expression(paste('Groundwater Flux  (', m^3, m^-2, d^-1, ')')),2,2.1,
          cex = 0.85)
    axis(1, at = seq(as.Date('2019-12-01'), as.Date('2020-03-01'), by = 'month'),
         labels = month.abb[c(12,1:3)], cex.axis = 0.8)
dev.off()



# gw correction?
d2 %>% select(site, date, group, DO.obs, ER) %>%
    left_join(select(ff, date, gw_flux_m3m2d)) %>%
    group_by(group) %>%
    summarize(across(c(DO.obs, gw_flux_m3m2d, ER), mean, na.rm = T)) %>%
    mutate(G = (DO.obs - 4) * gw_flux_m3m2d)
flux %>% filter(gas == 'CO2')%>%
ggplot(aes(-F_atm_gd, F_met_gd, col = delta_Q_md))+
    geom_point(size = 3)+
    geom_abline(slope = 1, intercept = 0)
flux %>% filter(gas == 'CO2')%>%
ggplot(aes(delta_Q_md, log(discharge), col = delta_Q_md))+
    geom_point(size = 3)+
    geom_abline(slope = 1, intercept = 0)

ggplot(flux, aes(distance_m, delta_Q_md, col = date, group = date)) +
  geom_line() +
  facet_wrap(~gas, scales = 'free_y', ncol = 1)

  # filter(site != "WB") %>%
ggplot(flux, aes(factor(date), F_atm_gm2d, fill = gas)) +
  geom_boxplot() +
  facet_wrap(~gas, scales = 'free_y', ncol = 1)+
  geom_hline(yintercept = 0)

flux_long <- flux %>%
  pivot_longer(cols = any_of(c('F_met_gm2d', 'F_atm_gm2d', 'F_missing_gm2d')),
               names_to = 'flux_source', values_to = 'flux_gm2d') %>%
  mutate(flux_source = factor(flux_source,
                              levels =c('F_met_gm2d','F_atm_gm2d',
                                                     'F_missing_gm2d') ))
# Figure 6 mass balance fluxies ####
png("figures/flux_components_by_site.png",
    width = 6, height = 6, res = 300, units = "in")
pal <- c("#009E73", "#56B4E9","#999999", "#E69F00" )
flux_long %>%
  # filter(flux_source != 'F_delta_gm2d')%>%
  mutate(flux_gm2d = case_when(reach_no == 2 & flux_source == 'F_missing_gm2d' ~
                                 NA_real_, TRUE ~flux_gm2d)) %>%
  group_by(reach_no, gas, flux_source) %>%
  summarize(flux_gm2d_sd = sd(flux_gm2d, na.rm = T),
            flux_gm2d_mean = mean(flux_gm2d, na.rm = T)) %>%
  # filter(reach_no != 2) %>%
ggplot(aes(reach_no, flux_gm2d_mean, fill = flux_source)) +
  geom_bar(stat = 'identity', position = 'dodge', width = .75) +
  geom_errorbar(aes(ymin = flux_gm2d_mean - flux_gm2d_sd,
                    ymax = flux_gm2d_mean + flux_gm2d_sd),
                position = position_dodge(.75),
                width = .2)+
  facet_wrap(~gas, scales = 'free_y', ncol = 1, strip.position = 'right')+
  geom_hline(yintercept = 0, lwd = .5) +
  ylab('Flux (g/m2/d)')+
  xlab('Stream reach (start - end, m)')+
  scale_fill_manual(values = pal,#c('forestgreen', 'lightblue', 'grey'),
                    name =  "Flux \nCategory",
                    labels = c("Aerobic \nMetabolism",
                                 "Atmospheric\nExchange ",
                                 "Missing \nSource/Sink")) +
  scale_x_discrete(limits = c(1:5),
                     # c('WBP','WB',
                     #          'CBP', 'PM','NHC'),
                   labels = c('0 - 2330', '2330 - 2500',
                              '2500 - 5000',
                              '5000 - 6880', '6880 - 8450'))+

  theme_bw()+
  theme(legend.text = element_text(margin = margin(t = 4,b = 4,  unit = "pt")))
dev.off()

png("figures/flux_components_over_time.png",
    width = 8, height = 6, res = 300, units = "in")
flux %>%
  as_tibble() %>%
  # mutate(F_met = F_met_gd/sa_m2,
  #        F_ev = F_ev_gd/sa_m2,
  #        F_miss = F_missing_gd/sa_m2) %>%
  pivot_longer(cols = any_of(c('F_met_gm2d', 'F_atm_gm2d', 'F_missing_gm2d')),
               names_to = 'flux_source', values_to = 'flux_gd') %>%
  mutate(sa_m2 = ifelse(!is.na(flux_gd), sa_m2, NA)) %>%
  group_by(date, gas, flux_source) %>%
  summarize(across(any_of(c('flux_gd', 'sa_m2')), ~sum(., na.rm = T)),
            flux_gm2d = flux_gd/sa_m2) %>%
  # filter(site != "WB") %>%
ggplot(aes(date, flux_gm2d, fill = flux_source)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gas, scales = 'free_y', ncol = 1)+
  geom_hline(yintercept = 0)
dev.off()


# water residence time vs gas residence time ####
res_time <- mb %>%
  mutate(vol_m3 = depth_m * sa_m2/length_m * 1000) %>%
  select(date, site, vol_m3) %>%
  left_join(flux, by = c('date', 'site')) %>%
  mutate(wrt_days = vol_m3/discharge /60 / 60/24,
         g_tot_g = vol_m3 * (c_in_ugl + c_out_ugl)/2/1000,
         grt_days = ifelse(gas == "O2", -g_tot_g/F_met_gd,
                           -g_tot_g/F_atm_gd)) %>%
  filter(!is.na(gas)) %>%
  select(date, site, reach_no, distance_m, discharge, gas, wrt_days, grt_days)
res_time %>%
  filter(site != 'WB', gas != 'O2') %>%
ggplot(aes(wrt_days, grt_days, col = date)) +
  geom_point(size = 2) +
  facet_grid(.~gas) +
  geom_abline(slope = 1, intercept = 0) +
  ylab('gas residence time (days)') +
  xlab('water residence time (days)')

write_csv(res_time, 'data/water_residence_times.csv')
flux %>%
  select(site, date, F_missing_gd, delta_Q_m3d, gas) %>%
  pivot_wider(names_from = gas, values_from = F_missing_gd) %>%
  ggplot(aes(date, CO2)) +
  geom_line() +
  facet_wrap(~site)
flux %>%
  mutate(delta = (Fout_gd - Fin_gd),
         metab = F_met_gd/delta,
         evasion = -F_atm_gd/delta,
         missing = F_missing_gd/delta,

         site = factor(site, levels = c('WBP','WB','CBP','PM','NHC'))) %>%
  # pivot_longer(cols = any_of(c('metab', 'evasion', 'missing')),
  #              names_to = 'flux_source', values_to = 'flux_gm2d') %>%
  # group_by(site, gas, flux_source) %>%
  # summarize(flux_sd = sd(flux_gm2d, na.rm = T),
  #           flux_gm2d = mean(flux_gm2d, na.rm = T)) %>%
  # ungroup() %>%
# filter(site != "WB") %>%
ggplot(aes(site, metab)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gas, scales = 'free_y', ncol = 1)+
  geom_hline(yintercept = 0)

flux %>%
  filter(site != "WB") %>%
ggplot(aes(delta_Q_m3d, F_missing_gd, col = c_out_ugl)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~gas, scales = 'free')

flux %>%
  filter(gas == "CO2") %>%
  mutate(CO2_gw_conc_mgl = ifelse(delta_Q_m3d >0 & F_missing_gd > 0,
                                  F_missing_gd/delta_Q_m3d, 0)) %>%
  ggplot(aes(site, CO2_gw_conc_mgl, col = site)) +geom_point()


# Old Calculations ####

dates <- unique(mb$date)
flux = data.frame()
# Calculate for each date:
for(i in 1:length(dates)){
  dd <- mb %>%
    filter(date == dates[i])
  n <- nrow(dd)-4
  tmp <- dd %>%
    slice(c(-(1:4))) %>%
    mutate(F_in_mgm2d = dd$discharge[1:n] *
             dd$concentration_ugl[1:n]/sa_m2*60*60*24,
           F_out_mgm2d = discharge * concentration_ugl/sa_m2*60*60*24,
           F_evas_mgm2d = (flux_ugld + dd$flux_ugld[1:n]) * depth_m / 2,
           F_instr_mgm2d = case_when(gas == "CO2" ~ (ER - GPP) * (44 / 32) *1000,
                                     gas == 'O2' ~ (GPP - ER) * 1000,
                                     TRUE ~ 0),
           F_ot_mgm2d =  F_out_mgm2d - F_in_mgm2d - F_instr_mgm2d +
             F_evas_mgm2d,
           gw_in_m3s = discharge - dd$discharge[1:n]) %>%
    select(c(site, date, gas, gw_in_m3s, ends_with("mgm2d")))

  dd <- dd %>% left_join(tmp, by = c("site", "date", "gas"))
  flux = bind_rows(flux, dd)
}

flux$site <- factor(flux$site, levels = c("UNHC", "WBP", "WB","CBP", "PM", "NHC"))

ggplot(flux, aes(x = site, F_ot_mgm2d)) +
  geom_boxplot() +
  facet_wrap(gas~., scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, col = "grey30")

longf <- flux %>%
  as.tibble() %>%
  mutate(F_evas_mgm2d = -F_evas_mgm2d) %>%
  pivot_longer(starts_with("F_"), names_to = "flux_category", values_to = "mgm2d")
png("figures/gas/flux_components_by_site.png",
    width = 8, height = 6, res = 300, units = "in")
  longf %>%
    filter(site != "UNHC",
           !(flux_category %in% c("F_in_mgm2d", "F_out_mgm2d"))) %>%
    ggplot(aes(site, (mgm2d), fill = flux_category)) +
      geom_boxplot(position = "dodge") +
      facet_wrap(gas~., scales = "free_y", ncol = 1) +
      geom_hline(yintercept = 0, col = "gray30")
dev.off()

ggplot(longf, aes(site, gw_in_m3s)) +#, col = factor(date))) +
  geom_boxplot()

longf %>%
  filter(site != "UNHC",
         flux_category == "F_ot_mgm2d") %>%
  ggplot(aes(gw_in_m3s/sa_m2 *60*60*24, mgm2d, col = site)) +
    geom_point() +
    facet_wrap(gas~., scales = "free_y") +
    xlab("ground water input (m/d)") +
    ylab("missing flux (mg/m2/d)") +
    geom_hline(yintercept = 0, col = "grey30")

write_csv(longf, "data/gas_data/mass_balance_nhc_ghg_data.csv")
chem <- dat %>%
  mutate(NH4NO3 = nh4n_mgl/no3n_mgl,
         CH4CO2 = (CH4.ugL/16)/(CO2.ugL/44),
         CH4O2 = (CH4.ugL/16)/(DO.obs/32)) %>%
  select(site, date, no3n_mgl, NH4NO3, CH4CO2, CH4O2, DO.obs, doc_mgl)

flux %>%
  left_join(chem, by = c("site", "date")) %>%
  ggplot(aes(CH4CO2, F_ot_mgm2d, col = site)) +
    geom_point() +
    facet_wrap(gas~., scales = "free_y", ncol = 1) +
    xlab("CH4:CO2") +
    xlim(0,.0087) +
    # xlim(0,1.8) +
    ylab("missing flux (mg/m2/d)") +
    geom_hline(yintercept = 0, col = "grey30")

png("figures/gas/missing_flux_by_gw_anaerobic_metric.png",
    width = 8, height = 7, res = 300, units = "in")
  flux %>%
    left_join(chem, by = c("site", "date")) %>%
    mutate('accumulated gw (m3/d)' = gw_in_m3s/60*60*24) %>%
    filter(NH4NO3 < 0.5,
           CH4CO2 < 0.0087) %>%
    pivot_longer(c("NH4NO3","CH4CO2", "accumulated gw (m3/d)"),
                 names_to = "element_ratio", values_to = "ratio_value") %>%
    ggplot(aes(ratio_value, (F_ot_mgm2d), col = site)) +
      geom_point() +
      facet_grid(gas~element_ratio, scales = "free") +
      xlab("") +
      # xlim(0,.5) +
      ylab("missing flux (mg/m2/d)") +
      geom_hline(yintercept = 0, col = "grey30")
dev.off()


ggplot(chem, aes(site, CH4O2)) +
  geom_boxplot()
-CO2_evas_mgm2d),col = "steelblue", lwd = 2) +
  geom_line(aes(y = CO2_gw_mgm2d), col = "sienna3", lwd = 2)# +
  geom_line(aes(y = CO2_ot_mgm2d), lty = 2, lwd = 2)# +
  legend()

# dc <- bind_rows(dd[1,], dc)
par(mfrow = c(3,1), mar = c(1,4,0,2), oma = c(4,0,4,0))
plot(dc$distance_m, dc$CO2_out_mgm2d - dc$CO2_in_mgm2d, lwd = 2,
     type = "l", ylim = c(-6000, 12000), ylab = "CO2 mg/m2/d",
      xaxt = 'n', col.bor = "grey30")
lines(dc$distance_m, dc$CO2_instr_mgm2d, col = "forestgreen", lwd = 2)
lines(dc$distance_m, dc$CO2_evas_mgm2d, col = "steelblue", lwd = 2)
lines(dc$distance_m, dc$CO2_gw_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
legend("topleft",
       legend = c("total flux", "In stream production", "evasion",
                  "missing flux (gw, anaerobic)"),
       lty = 1, lwd = 2, bty = 'n',
       col = c("black", "forestgreen", "steelblue","sienna3"))
mtext(dc$date[1], outer = T, cex = 1.2)

plot(dc$distance_m, dc$CH4_out_mgm2d - dc$CH4_in_mgm2d, lwd = 2,
     type = "l", ylim = c(-60, 30), ylab = "CH4 mg/m2/d",
     xaxt = 'n', col.bor = "grey30")
lines(dc$distance_m, dc$CH4_evas_mgm2d, col = "steelblue", lwd = 2)
lines(dc$distance_m, dc$CH4_ot_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
# legend("topleft",
#        legend = c("total flux", "evasion", "missing flux (gw, in stream)"),
#        lty = 1, lwd = 2, bty = 'n',
#        col = c("black", "steelblue","sienna3"))

plot(dc$distance_m, dc$N2O_out_mgm2d - dc$N2O_in_mgm2d, lwd = 2,
     type = "l", ylim = c(-.1, 1.5), ylab = "N2O mg/m2/d", xlab = "distance", bor.col = "grey30")
lines(dc$distance_m, dc$N2O_evas_mgm2d, col = "steelblue", lwd = 2)
lines(dc$distance_m, dc$N2O_ot_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
# legend("topleft",
#        legend = c("total flux", "evasion", "missing flux (gw, in stream)"),
#        lty = 1, lwd = 2, bty = 'n',
#        col = c("black", "steelblue","sienna3"))





ds <- data.frame(distance_m = seq(0, dc$distance_m[6], by = 10)) %>%
  left_join(dc, by = "distance_m") %>%
  mutate(across(c(-datetime, -date, -site), ~na.approx(.,na.rm = F)))

plot(ds$distance_m, ds$CO2_pred_ugL, type = "l", ylim = c(1, 7500), log = "y")

for(r in 2:nrow(ds)){
  lambda = (ds$K_CO2[r-1]/60/60/24)/
    (ds$discharge_m3s[r-1]/ds$depth_m[r-1]/ds$width_m[r-1])
  F_instr = (ds$ER[r] - ds$GPP[r])*(44/32) * ds$width_m[r] * 10
  F_out = ds$CO2_pred_ugL[r-1] * exp(-10* lambda) + F_instr/ds$discharge_m3s[r]/60/60/24
  ds$CO2_pred_ugL[r] <- F_out
}

lines(ds$distance_m, ds$CO2_pred_ugL, col = 2)

  g)
red_ugL = (CO2_in_mg + CO2_instr_mg - CO2_evas_mg)/
            discharge_m3s/60/60/24/ttime_d,
         CH4_pred_ugL = (CH4_in_mg - CH4_evas_mg)/
            discharge_m3s/60/60/24/ttime_d,
         N2O_pred_ugL = (N2O_in_mg - N2O_evas_mg)/
           discharge_m3s/60/60/24/ttime_d) %>%

ggplot(flux, aes(distance_m, CO2.ugL, color = factor(date))) +
  geom_point() + geom_line() + ylim(0,7500) +
  geom_point(aes(y = CO2_pred_ugL), color = 1)

plot(flux$distance_m, flux$CO2.ugL, type = "b")
points(flux$distance_m, flux$CO2_pred_ugL, col = "brown3")
t, distance_upstream_m, watertemp_C, GPP, ER,
         discharge_m3s, DO.obs, CH4.flux_ugld, CO2.flux_ugld, N2O.flux_ugld) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas", values_to = "flux_ugld") %>%
  mutate(gas = substr(gas, 1,3)) %>%
  left_join(gas, by = c("datetime", "site", "gas")) %>%
  mutate(date = as.Date(datetime),
         date = case_when(date == as.Date("2019-12-04") ~
                            as.Date("2019-12-03"),
                          date == as.Date("2020-01-30") ~
                            as.Date("2020-01-29"),
                          TRUE ~ date)) %>%
  arrange(date, distance_upstream_m)

flux$site <- factor(flux$site, levels = c("NHC", "PM", "CBP", "WB", "WBP","UNHC"))

# png("figures/gas/ghg_boxplots_bydate.png", height = 5, width = 7, units = "in", res = 300)
  ggplot(flux, aes(x = habitat, y = flux_ugld,
                   fill = habitat)) +
    geom_boxplot() +
    facet_wrap(gas~., scales = "free_y", ncol = 1)
# dev.off()
png("figures/gas/ghg_longitudinal_bydate.png", height = 5, width = 7, units = "in", res = 300)
ggplot(flux, aes(x = -distance_upstream_m, y = concentration_ugl,
                 color = factor(date))) +
  geom_line() +
  geom_point() +
  xlab("upstream --> downstream (m)")+
  facet_wrap(gas~., scales = "free_y", ncol = 1)
dev.off()
png("figures/gas/ghgflux_longitudinal_boxplots.png", height = 5, width = 4, units = "in", res = 300)
ggplot(flux, aes(x = site, y = flux_ugld))+
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(flux$site)))+
  xlab("upstream --> downstream (site)")+
  facet_wrap(gas~., scales = "free_y", ncol = 1)
dev.off()

ggplot(flux, aes(x = watertemp_C, concentration_ugl, color = date)) +
  geom_point(size = 2) +
  # geom_smooth(method = lm, color = "grey30") +
  # theme(legend.position = "bottom") +
  # scale_color_gradient(name = "log(Q)") +
  facet_wrap(gas~., scales = "free", ncol = 1)


# ccc <- ggplot(flux, aes(x = GPP+ER, y = concentration_ugl,
#                  color = watertemp_C)) +
#   geom_point(size = 2) +
#   geom_smooth(method = lm, color = "grey30") +
#   theme(legend.position = "bottom") +
#   scale_color_gradient(name = "temp C") +
#   facet_wrap(gas~., scales = "free_y", ncol = 1)
# fff <- ggplot(flux, aes(x = GPP+ER, y = flux_ugld,
#                  color = watertemp_C)) +
#   geom_point(size = 2) +
#   geom_smooth(method = lm, color = "grey30") +
#   theme(legend.position = "bottom") +
#   scale_color_gradient(name = "temp C") +
#   facet_wrap(gas~., scales = "free_y", ncol = 1)
# png("figures/gas/linear_ghg_by_NEPtemp.png", height = 7, width = 6, units = "in", res = 300)
# ggarrange(cc, ff)
# ggarrange(ccc, fff)
# dev.off()


# PCA ####
dat.pca <- d2 %>%
  select(-site, -habitat, -date, -K600, -ends_with("ugld")) %>%
  prcomp()

autoplot(dat.pca, data = d2, colour = "habitat")



# linear mixed effects models ####

# rescale covariates, normalize each to the mean
scaled <- dat %>%
  mutate(logQ = log(discharge_m3s),
         no3n_mgl = ifelse(no3n_mgl > 0.6, NA, no3n_mgl),
         doc_mgl = ifelse(doc_mgl < 0, 0, doc_mgl)) %>%
  select(site, habitat, logQ, watertemp_C, ER, GPP, DO.obs, no3n_mgl, doc_mgl,
         ends_with("ugld"), ends_with("ugL")) %>%
  mutate(across(-c(site, habitat), ~scale(.)[,1, drop = T]),
         across(c(site, habitat), ~factor(.)),
         ER = -ER,
         NEP = ER - GPP)
apply()
mm <- lmer(N2O.flux_ugld ~ CH4.flux_ugld + no3n_mgl + (1|site), data = scaled)

# steps_r <- mm_step$random %>%
#   as.tibble() %>%
#   mutate(predictor = rownames(mm_step$random),
#          effect = "random")
#
# steps_lme <- mm_step$fixed %>%
#   as.tibble() %>%
#   mutate(predictor = rownames(mm_step$fixed),
#          effect = "fixed") %>%
#   bind_rows(steps_r)

mm <- lmer(CH4.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs +
            CO2.ugL + (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER)
# r.squaredGLMM(mGPPER)
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

CH4_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "CH4")

mm <- lmer(CO2.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs +
             (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER)
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

CO2_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "CO2")

mm <- lmer(N2O.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs +
              (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER)
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

model_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "N2O") %>%
  bind_rows(CO2_steps,
            CH4_steps) %>%
  slice(c(-6, -7, -13, -14, -21, -22)) %>%
  select(gas, predictor, Eliminated, sum_sq = 'Sum Sq',
         NumDF, DenDF, F_value = 'F value', 'Pr(>F)')


write_csv(model_steps, "data/gas_data/Satterthwaite_DFmethod_lme_steps.csv")

# Analyze water chem data ####
spchem <- read_csv("data/water_chemistry/all_grab_data.csv") %>%
  filter(siteID %in% c("NHC", "UNHC")) %>%
  select(-flagID, -flagComment, -methodDetail, -writeInMethod, -regionID, -method) %>%
  group_by(siteID, dateTimeUTC, variable) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%#data.frame()
  pivot_wider(names_from = "variable", values_from = "value") %>%
  select(-TDN, -TN, -NH4, -PO4) %>%
  filter(across(-c(siteID, dateTimeUTC, TOC), ~!is.na(.)))

# pair with discharge and temperature
raw_dat <- data.frame()
for(site in c("NHC", "UNHC")){
  dat <- read_csv(paste0("data/metabolism/processed/", site, ".csv"),
                  guess_max = 10000)
  dat$site <- site
  raw_dat = bind_rows(raw_dat, dat)
}

chem <- raw_dat %>%
  as.tibble() %>%
  select(dateTimeUTC = DateTime_UTC, siteID = site, discharge, temp.water) %>%
  mutate(across(c(discharge, temp.water), na.approx, na.rm = T)) %>%
  right_join(spchem, by = c("dateTimeUTC", "siteID")) %>%
  filter(across(c(discharge, temp.water), ~!is.na(.))) %>%
  mutate(logQ = log(discharge),
         across(-c(siteID, dateTimeUTC, discharge), ~scale(.)[,1, drop = T]))

apply(chem, 2, function(x) sum(is.na(x)))

chem.pca <- chem %>%
  select(-siteID, -dateTimeUTC, -discharge, -TOC, -temp.water) %>%
  prcomp()

png("figures/gas/waterchem_pca.png", height = 5, width = 5, units = "in", res = 300)
autoplot(chem.pca, data = chem, colour = "siteID", size = 2,
         loadings = TRUE, loadings.label.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5) +
  # scale_color_gradient(low = "black", high = "red") +
  labs(title = "Water Chem by Discharge @ NHC, UNHC")
dev.off()
# biplot(chem.pca, pch = 20)
