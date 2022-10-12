# Plot instantaneous sensor data used to calculate metabolism:
library(tidyverse)
library(lubridate)
library(ggpubr)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")


dat <- read_csv('data/site_data/metabolism.csv')


d <- dat %>%
    mutate(sitename = case_when(site == 'NHC' ~ 'NHC_8.5',
                                site == 'PM' ~ 'NHC_6.9',
                                site == 'CBP' ~ 'NHC_5',
                                site == 'WB' ~ 'NHC_2.5',
                                site == 'WBP' ~ 'NHC_2.3',
                                site == 'UNHC' ~ 'NHC_0'),
           site = factor(sitename, levels = c("NHC_0", "NHC_2.3", "NHC_2.5",
                                              "NHC_5", "NHC_6.9", "NHC_8.5")),
           datetime = with_tz(DateTime_UTC, 'EST')) %>%
    select(datetime, site, discharge, temp.water, DO.obs)


DO <- ggplot(d, aes(datetime, DO.obs, group = site, col = site))+
    geom_line() +
    theme_bw() +
    ylab(expression(paste('Dissolved ', O[2], ' (mg ',l^-1,')' ))) +
    xlab('Date')
# DO <- DO +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           plot.margin = unit(c(0.5,0,0,0), 'lines'))
wtemp <- ggplot(d, aes(datetime, temp.water, group = site, col = site))+
    geom_line() +
    theme_bw()+
    ylab(expression(paste('Water Temp (',degree,'C)' ))) +
    xlab('Date')
# wtemp <- wtemp +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           plot.margin = unit(c(0.2,0,0,0), 'lines'))
discharge <- ggplot(d, aes(datetime, discharge, group = site, col = site))+
    geom_line() +
    theme_bw() +
    scale_y_log10()+
    ylab(expression(paste('Discharge (',m^3, s^-1,')' ))) +
    xlab('Date')
# discharge <- discharge +
#     theme(plot.margin = unit(c(0.5,0.1,0.2,0), 'lines'))

tiff('figures/final/raw_metab_inputs.tiff', width = 9, height = 5.5,
     res = 300, units = 'in')
ggpubr::ggarrange(DO + rremove('xlab') + rremove('x.text'),
                  wtemp + rremove('xlab') + rremove('x.text'), discharge,
                  ncol = 1, common.legend = TRUE, legend = 'right', align = 'hv')
dev.off()
