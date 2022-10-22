# plot fraction of CO2 produced instream
setwd('C:/Users/alice.carter/git/ghg_patterns_nhc/')
# devtools::install_github("psyteachr/introdataviz")
library(tidyverse)
library(ggpubr)
library(introdataviz)
library(viridis)

dat <- read_csv('data/fraction_of_instream_production_CO2_and_CH4CO2ratios.csv') %>%
    mutate(distance_m = 8450 - distance_upstream_m,
           distance_m = factor(distance_m))
dvs <- read_csv("data/ghg_filled_drivers_dataframe.csv") %>%
    filter(site == 'NHC')
           # date >= min(unique(dat$group)),
           # date <= max(unique(dat$group)))

gw <- read_csv('data/gw_fluxes.csv')

dat%>%
    select(site, date, group, CO2_flux, NEP_CO2, instr, ER, GPP) %>%
    left_join(gw) %>% select(-date) %>%
    # filter(CO2_flux >0) %>%
    mutate(gw_corr = case_when(gw_flux_m3m2d > 0 ~ gw_flux_m3m2d * 8,
                               TRUE ~ 0),
           ER_corr = ER - gw_corr,
           ER_corr = ifelse(ER_corr < 0, 0, ER_corr),
           NEP_corr = ER_corr - GPP,
           NEP_CO2_corr = ER_corr * 44/32,
           pctNEP = NEP_CO2/CO2_flux,
           pctNEP_corr = NEP_CO2_corr/CO2_flux,
           pctNEP_filt = case_when(pctNEP < 0 ~ 0,
                                   pctNEP >1 ~1,
                                   TRUE ~pctNEP)) %>% #filter(pctNEP_filt == 0)
    summary()
dat%>%
    select(site, group, CO2_flux, NEP_CO2, instr) %>%
    mutate(pctNEP = NEP_CO2/CO2_flux,
           pctNEP_filt = case_when(pctNEP < 0 ~ 0,
                                   pctNEP >1 ~1,
                                   TRUE ~pctNEP)) %>% #filter(pctNEP_filt == 0)
    group_by(site) %>%# filter(CO2_flux >0) %>%
    summarize(mean= mean(pctNEP_filt, na.rm = T),
              median= median(pctNEP_filt, na.rm = T),
              sd = sd(pctNEP_filt, na.rm = T),
              min = min(pctNEP_filt, na.rm = T) ,
              xs = max(pctNEP, na.rm = T),
              all = length(which(pctNEP_filt == 1)),
              zero = length(which(pctNEP_filt == 0)))
dat%>%
    select(site, group, CO2_flux, NEP_CO2, instr) %>%
    # filter(CO2_flux >0) %>%
    mutate(pctNEP = NEP_CO2/CO2_flux,
           pctNEP_filt = case_when(pctNEP < 0 ~ 0,
                                   pctNEP >1 ~1,
                                   TRUE ~pctNEP)) %>% #filter(pctNEP_filt == 0)
    summarize(all = sum(!is.na(CO2_flux)),
              one = length(which(pctNEP_filt == 1)),
              zero = length(which(pctNEP_filt == 0)),
              a = length(which(NEP_CO2<0 &CO2_flux >0)),
              b = length(which(NEP_CO2>0 &CO2_flux <0)),
              p_one = one/all,
              p_zero = zero/all)
dat%>%
    select(site, group, CO2_flux, NEP_CO2, instr) %>%
    mutate(pctNEP = NEP_CO2/CO2_flux,
           pctNEP_filt = case_when(pctNEP < 0 ~ 0,
                                   pctNEP >1 ~1,
                                   TRUE ~pctNEP)) %>%
    mutate(fall = ifelse(group > ymd('2020-01-15'), 'spring', 'fall')) %>%
    group_by(fall) %>%
    summarize(mean= mean(pctNEP_filt, na.rm = T),
              sd = sd(pctNEP_filt, na.rm = T),
              min = min(pctNEP_filt, na.rm = T) ,
              xs = max(pctNEP, na.rm = T),
              all = length(which(pctNEP_filt == 1)),
              zero = length(which(pctNEP_filt == 0)))
    ggplot(aes(group, pctNEP_filt, col = site)) +
    geom_point()

dat %>%
    select(CO2_flux, NEP_CO2) %>%
    filter(CO2_flux >0) %>%
    mutate(NEP_CO2 = ifelse(NEP_CO2>0, NEP_CO2, 0)) %>%
    summarize(across(everything(), sum))

dat%>%
    group_by(group)%>%
    summarize(CO2_flux = median(CO2_flux, na.rm = T)) %>%
    ggplot(aes(group, CO2_flux)) +
    geom_point()
plot(dat$date, dat$CO2_flux)


d1 <- dat %>%
    filter(CO2_flux > 0)%>%
    select(distance_m, date, CO2_flux, NEP_CO2) %>%
    pivot_longer(cols = any_of(c("NEP_CO2", "CO2_flux")),
                 values_to = "gm2d",names_to = "category")
d2 <- dat %>%
    filter(CO2_flux > 0)%>%
    select(distance_m, date, CO2_flux, NEP_CO2 = NEP_CO2_low) %>%
    pivot_longer(cols = any_of(c("NEP_CO2", "CO2_flux")),
                 values_to = "gm2d_low",names_to = "category") %>%
    mutate(gm2d_low = ifelse(category == 'NEP_CO2', gm2d_low, NA_real_))
d3 <- dat %>%
    filter(CO2_flux > 0)%>%
    select(distance_m, date, CO2_flux, NEP_CO2 = NEP_CO2_high) %>%
    pivot_longer(cols = any_of(c("NEP_CO2", "CO2_flux")),
                 values_to = "gm2d_high",names_to = "category") %>%
    mutate(gm2d_high = ifelse(category == 'NEP_CO2', gm2d_high, NA_real_))
dd <- d1 %>% full_join(d2) %>% full_join(d3) %>%
    mutate(across(starts_with('gm2d'), ~ifelse(.>0, ., NA_real_)))
dmeans <- dd %>% group_by(distance_m, category)%>%
    summarize(across(starts_with('gm2'), median, na.rm = T)) %>%
    arrange(category)
dmeans$gm2d_high[1:6] <-
    dmeans$gm2d_low[1:6] <-
    dmeans$gm2d[1:6]

col.NEP <- 'forestgreen'
col.NER <- plasma(7)[5]
col.xs <- 'grey80'
dp <- dat %>%
    group_by(group) %>%
    summarize(CO2_flux = median(CO2_flux, na.rm = T),
              NEP_CO2 = median(NEP_CO2, na.rm = T)) %>%
    mutate(CO2_neg = ifelse(CO2_flux < 0, CO2_flux, 0),
           NEP_neg = case_when(NEP_CO2 < 0 & NEP_CO2 > CO2_flux ~ NEP_CO2,
                               NEP_CO2 >= 0 ~0,
                               NEP_CO2 < CO2_flux ~ 0),
           CO2_pos = ifelse(CO2_flux>0, CO2_flux, 0),
           NEP_pos = case_when(NEP_CO2 > 0 & NEP_CO2 < CO2_flux ~ NEP_CO2,
                               NEP_CO2 <= 0 ~0,
                               TRUE ~ CO2_flux))
png('figures/CO2_flux_from_NEP_date_polygon.png', width = 6.5, height = 4,
    res = 300,  units = 'in')
ylim = range(c(dp$CO2_flux, dp$NEP_CO2)) *1.1
plot(dp$group, dp$CO2_flux, type = 'n',
     ylab = expression(paste(CO[2], ' flux (g/', m^2, '/d)')),
     xlab = 'Sample Date', ylim = ylim)
polygon(c(dp$group, rev(dp$group)), c(dp$CO2_pos, rev(dp$NEP_pos)),
        col = col.xs, border = NA)
polygon(c(dp$group, rev(dp$group)), c(rep(0, nrow(dp)), rev(dp$NEP_pos)),
        col = alpha(col.NER, 0.8), border = NA)
polygon(c(dp$group, rev(dp$group)), c(dp$CO2_neg, rev(dp$NEP_neg)),
        col = col.xs, border = NA)
polygon(c(dp$group, rev(dp$group)), c(rep(0, nrow(dp)), rev(dp$NEP_neg)),
        col = alpha(col.NER, .8), border = NA)
abline(v = unique(dp$group), lty = 2)
polygon(c(dp$group, rev(dp$group)), c(dp$CO2_pos, rep(4.85, nrow(dp))),
        col = 'white', border = NA)
polygon(c(dp$group, rev(dp$group)), c(dp$CO2_neg, rep(-0.45, nrow(dp))),
        col = 'white', border = NA)
dev.off()

# png('figures/CO2_flux_from_NEP_site_violins.png', width = 6.5, height = 4,
#     res = 300,  units = 'in')
viols <- ggplot(dd, aes(distance_m, gm2d, fill = category) )+
    introdataviz::geom_split_violin(alpha = .4, trim = FALSE, adjust = 1.1,
                                    width = 1.4, lwd = .25) +
    introdataviz::geom_split_violin(aes(y = gm2d_low),
                                    alpha = .4, trim = FALSE, adjust = 1.1,
                                    width = 1.4, lwd = .25) +
    introdataviz::geom_split_violin(aes(y = gm2d_high),
                                    alpha = .4, trim = FALSE, adjust = 1.1,
                                    width = 1.4, lwd = .25) +
    geom_boxplot(data = dmeans, aes(x = distance_m, y = gm2d,
                             lower = gm2d, upper = gm2d, middle = gm2d,
                             fill = category),
                 width = .2, alpha = .6, fatten = NULL, show.legend = FALSE,
                 outlier.shape = NA, coef = 0) +
    geom_boxplot(data = dmeans, aes(x = distance_m, y = gm2d_low,
                             lower = gm2d_low, upper = gm2d_low, middle = gm2d_low,
                             fill = category),
                 width = .2, alpha = .6, fatten = NULL, show.legend = FALSE,
                 outlier.shape = NA, coef = 0) +
    geom_boxplot(data = dmeans, aes(x = distance_m, y = gm2d_high,
                             lower = gm2d_high, upper = gm2d_high, middle = gm2d_high,
                             fill = category),
                 width = .2, alpha = .6, fatten = NULL, show.legend = FALSE,
                 outlier.shape = NA, coef = 0) +
    scale_fill_manual(values = c(col.xs, col.NER),
                      name = '',
                      labels = c(expression(paste('Total ', CO[2], ' flux  ')),
                                 expression(paste(CO[2], ' from NEP')),
                                 'median'))+
    labs(x = 'Distance Downstream (m)',
         y = expression(paste(CO[2], ' flux (g'~m^-2~d^-1*')')))+
    theme_bw()+
    scale_x_discrete(expand = c(.18,0)) +
    theme(legend.position = 'top',
          axis.title = element_text(size = 9),
          legend.key.size = unit(0.3,'cm'),
          legend.text = element_text(size = 9),
          # legend.direction = 'vertical',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

viols_1d <- ggplot(dd, aes(distance_m, gm2d, fill = category) )+
    introdataviz::geom_split_violin(alpha = .8, trim = FALSE, adjust = 1.1,
                                    width = 1.4, lwd = .25) +
    geom_boxplot(data = dmeans, aes(x = distance_m, y = gm2d,
                             lower = gm2d, upper = gm2d, middle = gm2d,
                             fill = category),
                 width = .3, alpha = .6, fatten = NULL, show.legend = FALSE,
                 outlier.shape = NA, coef = 0) +
    scale_fill_manual(values = c(col.xs, col.NER),
                      name = '',
                      labels = c(expression(paste('Total ', CO[2], ' flux  ')),
                                 expression(paste(CO[2], ' from NEP')),
                                 'median'))+
    labs(x = 'Distance Downstream (m)',
         y = expression(paste(CO[2], ' flux (g'~m^2~d^-1*')')))+
    theme_bw()+
    scale_x_discrete(expand = c(.18,0)) +
    theme(legend.position = 'top',
          axis.title = element_text(size = 9),
          legend.key.size = unit(0.3,'cm'),
          legend.text = element_text(size = 9),
          # legend.direction = 'vertical',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
# dev.off()
# png('figures/CO2_flux_from_NEP_date_bar.png', width = 6.5, height = 4,
#     res = 300,  units = 'in')
dpdevs <- dp %>%
    mutate(NEP_CO2_high = NEP_pos/0.8,
           NEP_CO2_low = NEP_pos/0.8 * 0.6,
           NEP_CO2_high = ifelse(NEP_CO2_high > CO2_pos, CO2_pos, NEP_CO2_high))

bars_var <- dpdevs %>%
    mutate(NEP_pos = NEP_pos + NEP_neg,
        Extra = CO2_flux - NEP_pos)%>%
    pivot_longer(cols = any_of(c('NEP_pos', 'Extra')),
                 names_to = 'category', values_to = 'gm2d') %>%
    ggplot(aes(group, gm2d, fill = category) )+
    geom_bar(stat = 'identity') + #ylim(-.1,1) +
    geom_errorbar(aes(x = group, ymax = NEP_CO2_high,
                                     ymin = NEP_CO2_low),
                  position = 'identity', width = 2) +
    geom_hline(yintercept = 0, size = .3)+
    xlab('Sample Date')+
    guides(fill = F)+
    scale_fill_manual(values = c(col.xs, alpha(col.NER, .8)))+
    scale_y_continuous(name = expression(paste(CO[2], ' flux (g'~m^2~d^-1*')')),
                       sec.axis = sec_axis(trans = ~.*3,
                                name = expression(paste('Discharge (',
                                                        m^3~s^-1*')'))))+
    geom_line(data = dvs,
              aes(x = date, y = ((discharge))/3, fill = NULL),
              size = .75, col = 'steelblue4')+
    theme_bw()+
    theme(plot.margin = unit(c(.1,.2,.2,.2), "cm"),
          axis.title = element_text(size = 9),
          axis.title.y.right = element_text(color = 'steelblue4', size = 9.5),
          axis.text.y.right = element_text(color = 'steelblue4'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
bars <- dp %>%
    mutate(NEP_pos = NEP_pos + NEP_neg,
        Extra = CO2_flux - NEP_pos)%>%
    pivot_longer(cols = any_of(c('NEP_pos', 'Extra')),
                 names_to = 'category', values_to = 'gm2d') %>%
    ggplot(aes(group, gm2d, fill = category) )+
    geom_bar(stat = 'identity') + #ylim(-.1,1) +
    geom_hline(yintercept = 0, size = .3)+
    xlab('Sample Date')+
    guides(fill = F)+
    scale_fill_manual(values = c(col.xs, alpha(col.NER, .8)))+
    scale_y_continuous(name = expression(paste(CO[2], ' flux (g'~m^2~d^-1*')')),
                       sec.axis = sec_axis(trans = ~.*3,
                                name = expression(paste('Discharge (',
                                                        m^3~s^-1*')'))))+
    geom_line(data = dvs,
              aes(x = date, y = ((discharge))/3, fill = NULL),
              size = .75, col = 'steelblue4')+
    theme_bw()+
    theme(plot.margin = unit(c(.1,.2,.2,.2), "cm"),
          axis.title = element_text(size = 9),
          axis.title.y.right = element_text(color = 'steelblue4', size = 9.5),
          axis.text.y.right = element_text(color = 'steelblue4'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
# dev.off()

tiff('figures/final/CO2_flux_from_NEP_combined.tif', width = 6.5, height = 3,
    res = 300, units = 'in')
    ggarrange(viols, bars_var, labels = c('a', 'b'), align = 'h',
              legend = 'top', common.legend = TRUE, widths = c(1,1.3))
dev.off()
tiff('figures/final/CO2_flux_from_NEP_combined_1d.tif', width = 6.5, height = 3,
    res = 300, units = 'in')
    ggarrange(viols_1d, bars, labels = c('a', 'b'), align = 'h',
              legend = 'top', common.legend = TRUE, widths = c(1,1.3))
dev.off()
