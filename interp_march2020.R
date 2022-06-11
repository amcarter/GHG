library(tidyverse)
library(lubridate)
library(mice)

setwd('C:/Users/alice.carter/git/nhc_50yl/')
# source('src/helpers.R')
filelist <- list.files('data/metabolism/processed/')
dat <- tibble()
for(i in filelist){
    dd <- read_csv(paste0('data/metabolism/processed/', i))
    dd$site <- dd$site[1]
    dat <- bind_rows(dat, dd)
}
dat <- filter(dat, site!= 'PWC')
dat <- filter(dat, !(site == 'PM' & DateTime_UTC > ymd_hms('2020-03-07 05:00:00')),
              DateTime_UTC < ymd_hms('2020-03-24 05:00:00'))
# Oxygen ####
d = dat %>%
    select(DateTime_UTC, site, DO.obs, discharge) %>%
    pivot_wider(names_from = site,
                values_from = -all_of(c('DateTime_UTC', 'site'))) %>%
    select(DateTime_UTC, starts_with('DO.obs'))

#initial inspection of one gap section
d_plotinds = which(d$DateTime_UTC > as.POSIXct('2020-02-20', tz = 'UTC') &
                       d$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
par(mfrow=c(3, 2), mar=c(1,0,1,0), oma=c(0,0,0,0))
for(dfc in colnames(d)[-1]){
    # plot(d$DateTime_UTC, d[, dfc, drop = TRUE],
    plot(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='orange', type='p', lwd=2, pch = '.')
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-2)
}

#seasonally splitted imputation ####
d_interp <- d %>%
    mutate(across(-DateTime_UTC, ~imputeTS::na_seasplit(ts(., deltat = 1/96))))

par(mfrow=c(3, 2), mar=c(1,0,1,0), oma=c(0,0,0,0))
for(dfc in colnames(d_interp)[-1]){
    plot(d_interp$DateTime_UTC[d_plotinds], d_interp[d_plotinds, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    lines(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
          col='black', lwd=1)
    mtext(dfc, 3, line=-2)
}

d_sub <- filter(d, DateTime_UTC >= ymd_hms('2020-02-20 05:00:00'))
dd <- d_interp %>%
    # select(DateTime_UTC, DO.obs_WB, DO.obs_WBP, DO.obs_NHC) %>%
    filter(DateTime_UTC >= ymd_hms("2020-02-20 05:00:00"))

nhc <- ts(dd$DO.obs_NHC, frequency = 96, start = d_sub$DateTime_UTC[1])
nhcd <- decompose(nhc)

unhc <- ts(dd$DO.obs_UNHC, frequency = 96, start = d_sub$DateTime_UTC[1])
unhcd <- decompose(unhc)

# plot(nhcd$trend)
# lines(unhcd$trend, lty = 2)

w <- which(is.na(d_sub$DO.obs_UNHC))
w <- c(w[1]-1, w, w[length(w)] + 1)
a <- unhcd$trend[w[1]] - nhcd$trend[w[1]]
b <- unhcd$trend[w[length(w)]] - nhcd$trend[w[length(w)]]
line <- seq(a, b, length.out = length(w))
unhct <- unhcd$trend
unhct[w] <- nhcd$trend[w] + line
# lines(unhct, lty = 2, col = 3)
dd <- dd %>%
    mutate(DO.obs_UNHC = DO.obs_UNHC - unhcd$trend + unhct)

pm <- ts(dd$DO.obs_PM, frequency = 96, start = d_sub$DateTime_UTC[1])
pm <- decompose(pm)

w <- which(is.na(d_sub$DO.obs_PM))
w <- c(w[1]-1, w)
a <- pm$trend[w[1]] - nhcd$trend[w[1]]
pmt <- pm$trend
pmt[w] <- nhcd$trend[w] + a

dd <- dd %>%
    mutate(DO.obs_PM = DO.obs_PM - pm$trend + pmt)

dd$DO.obs_CBP[is.na(d_sub$DO.obs_CBP)] <- NA
par(mfrow=c(3, 2), mar=c(1,0,1,0), oma=c(0,0,0,0))
for(dfc in colnames(dd)[-1]){
    plot(dd$DateTime_UTC, dd[, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    lines(d_sub$DateTime_UTC, d_sub[, dfc, drop = TRUE],
          col='black', lwd=1)
    mtext(dfc, 3, line=-2)
}

dd_DO <- dd

# Temperature ####
d = dat %>%
    select(DateTime_UTC, site, temp.water, DO.obs) %>%
    pivot_wider(names_from = site,
                values_from = -all_of(c('DateTime_UTC', 'site'))) %>%
    select(DateTime_UTC, starts_with('temp.water'))

#initial inspection of one gap section
d_plotinds = which(d$DateTime_UTC > as.POSIXct('2020-02-20', tz = 'UTC') &
                       d$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
par(mfrow=c(3, 2), mar=c(1,0,1,0), oma=c(0,0,0,0))
for(dfc in colnames(d)[-1]){
    # plot(d$DateTime_UTC, d[, dfc, drop = TRUE],
    plot(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='orange', type='p', lwd=2, pch = '.')
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-2)
}

lm(temp.water_PM ~ temp.water_NHC + temp.water_UNHC, data = d)
w <- which(is.na(d$temp.water_PM))
d$temp.water_PM[w] <- -1.3777 + 0.2658 * d$temp.water_NHC[w] +
    0.8904 * d$temp.water_UNHC[w]

lm(temp.water_WBP ~ temp.water_NHC + temp.water_UNHC, data = d)
w <- which(is.na(d$temp.water_WBP))
d$temp.water_WBP[w] <- -.2599 + 0.587 * d$temp.water_NHC[w] +
    0.4566 * d$temp.water_UNHC[w]

lm(temp.water_WB ~ temp.water_NHC + temp.water_UNHC, data = d)
w <- which(is.na(d$temp.water_WB))
d$temp.water_WB[w] <- -.4442 + 0.2514 * d$temp.water_NHC[w] +
    0.794 * d$temp.water_UNHC[w]

dd <- filter(d, DateTime_UTC >= ymd_hms("2020-02-20 05:00:00"))

par(mfrow=c(3, 2), mar=c(1,0,1,0), oma=c(0,0,0,0))
for(dfc in colnames(dd)[-1]){
    plot(dd$DateTime_UTC, dd[, dfc, drop = TRUE],
         col='orange', type='l', lwd=2)
    mtext(dfc, 3, line=-2)
}

airpres <- read_csv('data/siteData/NOAA_airpres.csv')
dd_temp <- dd
dd <- full_join(dd_DO, dd_temp)
dd <- left_join(dd, airpres) %>%
    mutate(DO.sat_PM = streamMetabolizer::calc_DO_sat(temp.water_PM,
                                                      air_kPa * 10),
           DO.sat_UNHC = streamMetabolizer::calc_DO_sat(temp.water_UNHC,
                                                      air_kPa * 10),
           DO.sat_NHC = streamMetabolizer::calc_DO_sat(temp.water_NHC,
                                                      air_kPa * 10),
           DO.sat_WB = streamMetabolizer::calc_DO_sat(temp.water_WB,
                                                      air_kPa * 10),
           DO.sat_WBP = streamMetabolizer::calc_DO_sat(temp.water_WBP,
                                                      air_kPa * 10))

write_csv(dd, '../ghg_patterns_nhc/data/processed_sensor/interp_DO_march2020.csv')

# replace in processed files:
# PM:
d <- dd %>% select(DateTime_UTC, ends_with('PM'))
dat <- read_csv("data/metabolism/processed/PM.csv", guess_max = 10000)

dat <- full_join(dat, d, by = 'DateTime_UTC')
dat <- dat %>%
    mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs_PM),
           temp.water = ifelse(!is.na(temp.water),
                                temp.water, temp.water_PM),
           DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat_PM)) %>%
    select(-ends_with('_PM'))
write_csv(dat, "data/metabolism/processed/PM.csv")
# UNHC:
d <- dd %>% select(DateTime_UTC, ends_with('UNHC'))
dat <- read_csv("data/metabolism/processed/UNHC.csv", guess_max = 10000)

dat <- full_join(dat, d, by = 'DateTime_UTC')
dat <- dat %>%
    mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs_UNHC),
           temp.water = ifelse(!is.na(temp.water),
                                temp.water, temp.water_UNHC),
           DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat_UNHC)) %>%
    select(-ends_with('_UNHC'))
write_csv(dat, "data/metabolism/processed/UNHC.csv")
# NHC:
d <- dd %>% select(DateTime_UTC, ends_with('NHC'))
dat <- read_csv("data/metabolism/processed/NHC.csv", guess_max = 10000)

dat <- full_join(dat, d, by = 'DateTime_UTC')
dat <- dat %>%
    mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs_NHC),
           temp.water = ifelse(!is.na(temp.water),
                                temp.water, temp.water_NHC),
           DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat_NHC)) %>%
    select(-ends_with('_NHC'))
write_csv(dat, "data/metabolism/processed/NHC.csv")
# WBP:
d <- dd %>% select(DateTime_UTC, ends_with('WBP'))
dat <- read_csv("data/metabolism/processed/WBP.csv", guess_max = 10000)

dat <- full_join(dat, d, by = 'DateTime_UTC')
dat <- dat %>%
    mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs_WBP),
           temp.water = ifelse(!is.na(temp.water),
                                temp.water, temp.water_WBP),
           DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat_WBP)) %>%
    select(-ends_with('_WBP'))
write_csv(dat, "data/metabolism/processed/WBP.csv")
# WB:
d <- dd %>% select(DateTime_UTC, ends_with('WB'))
dat <- read_csv("data/metabolism/processed/WB.csv", guess_max = 10000)

dat <- full_join(dat, d, by = 'DateTime_UTC')
dat <- dat %>%
    mutate(DO.obs = ifelse(!is.na(DO.obs), DO.obs, DO.obs_WB),
           temp.water = ifelse(!is.na(temp.water),
                                temp.water, temp.water_WB),
           DO.sat = ifelse(!is.na(DO.sat), DO.sat, DO.sat_WB)) %>%
    select(-ends_with('_WB'))
write_csv(dat, "data/metabolism/processed/WB.csv")

