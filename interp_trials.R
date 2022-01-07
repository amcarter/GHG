library(tidyverse)
library(lubridate)
library(mice)

setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/ghg_patterns_nhc/')
# source('src/helpers.R')

d = read_csv('data/processed_sensor/compiled_nhc_dat.csv')

d = d %>%
    select(DateTime_UTC, site, DO.obs, discharge) %>%
    pivot_wider(names_from = site,
                values_from = -all_of(c('DateTime_UTC', 'site'))) %>%
    select(DateTime_UTC, starts_with('DO.obs'), discharge_NHC)

#initial inspection of one gap section
d_plotinds = which(d$DateTime_UTC > as.POSIXct('2020-03-01', tz = 'UTC') &
                       d$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in colnames(d)[-1]){
    # plot(d$DateTime_UTC, d[, dfc, drop = TRUE],
    plot(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='orange', type='p', lwd=2, pch = '.')
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-4)
}

#seasonally splitted imputation ####
d_interp <- d %>%
    mutate(across(-DateTime_UTC, ~imputeTS::na_seasplit(ts(., deltat = 1/96))))

par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in colnames(d_interp)[-1]){
    plot(d_interp$DateTime_UTC[d_plotinds], d_interp[d_plotinds, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    lines(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
          col='black', lwd=1)
    mtext(dfc, 3, line=-4)
}

dd <- d_interp %>% 
    select(DateTime_UTC, DO.obs_WB, DO.obs_WBP, DO.obs_UNHC) %>% 
    filter(DateTime_UTC > ymd_hms("2020-03-01 00:00:00"))
write_csv(dd, 'data/processed_sensor/interp_DO_march2020.csv')

# replace in processed files:
    site = 'WBP'
    dat <- read_csv(paste0("../hall_50yl2/NHC_2019_metabolism/data/metabolism/processed/", 
                           site, ".csv"), guess_max = 10000)
    dat <- full_join(dat, tmp, by = 'DateTime_UTC') %>%
        mutate(DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat2),
               temp.water = ifelse(!is.na(temp.water), 
                                    temp.water, watertemp_C2)) %>%
        select(-DO.sat2, -watertemp_C2)
    write_csv(dat, paste0("../hall_50yl2/NHC_2019_metabolism/data/metabolism/processed/", 
                          site, ".csv"))

# make sure they have all the other business they need:
tmp = read_csv('data/processed_sensor/compiled_nhc_dat.csv') %>%
    filter(DateTime_UTC > ymd_hms("2020-03-01 00:00:00")) %>%
    group_by(DateTime_UTC) %>%
    summarize(DO.sat2 = mean(DO.sat, na.rm = T), 
              watertemp_C2 = mean(watertemp_C, na.rm = T))

ggplot(dat, aes(DateTime_UTC, watertemp_C)) +
    geom_line()

    
# mice (one of 5 models) ####

d_mice <- mice(d[, -1], m=5, method = 'pmm')
d_list <- mice::complete(d_mice, action = 'all')
d_imp_mat <- lapply(d_list, as.matrix)
d_ave <- (d_imp_mat[[1]] + d_imp_mat[[2]] + d_imp_mat[[3]] + d_imp_mat[[4]] + d_imp_mat[[5]]) / 5
d_imp <- bind_cols(select(d, DateTime_UTC), as_tibble(d_ave))

# model_fit <- with(data = d_mice, exp = lm(DO.obs_NHC ~ DO.obs_UNHC))
# d_mice_pooled <- pool(model_fit)

par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in colnames(d_imp)[-1]){
    # plot(d$DateTime_UTC, d[, dfc, drop = TRUE],
    plot(d_imp$DateTime_UTC[d_plotinds], d_imp[d_plotinds, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    lines(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='black', lwd=1)
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-4)
}

# mice (ensemble approach) ####

d_1 <- mice::complete(d_mice, action = 1) %>%
    bind_cols(select(d, DateTime_UTC))
d_2 <- mice::complete(d_mice, action = 2) %>%
    bind_cols(select(d, DateTime_UTC))
d_3 <- mice::complete(d_mice, action = 3) %>%
    bind_cols(select(d, DateTime_UTC))
d_4 <- mice::complete(d_mice, action = 4) %>%
    bind_cols(select(d, DateTime_UTC))
d_5 <- mice::complete(d_mice, action = 5) %>%
    bind_cols(select(d, DateTime_UTC))

par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in colnames(d_1)[-1]){
    plot(d_1$DateTime_UTC[d_plotinds], d_1[d_plotinds, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    lines(d_2$DateTime_UTC[d_plotinds], d_2[d_plotinds, dfc, drop = TRUE],
         col='red', lwd=1)
    lines(d_3$DateTime_UTC[d_plotinds], d_3[d_plotinds, dfc, drop = TRUE],
         col='purple', lwd=1)
    lines(d_4$DateTime_UTC[d_plotinds], d_4[d_plotinds, dfc, drop = TRUE],
         col='blue', lwd=1)
    lines(d_5$DateTime_UTC[d_plotinds], d_5[d_plotinds, dfc, drop = TRUE],
         col='darkgreen', lwd=1)
    lines(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='black', lwd=1)
    mtext(dfc, 3, line=-4)
}
