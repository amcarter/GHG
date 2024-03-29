# Initial GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(RColorBrewer)
library(tidyverse)
library(lubridate)
library(ggfortify)
library(ggpubr)
library(agricolae)
library(zoo)
setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")

dat <- read_csv("data/ghg_flux_complete_drivers_dataframe.csv")%>%
  filter(!is.na(datetime),
         site != "MC751")
dvs <- read_csv("data/ghg_filled_drivers_dataframe.csv")

ggplot(dat, aes(site, K_CO2, col = factor(site))) +
  geom_boxplot()
d2 <- read_csv("data/ghg_flux_complete_drivers_dataframe_noNAs.csv")
# plots ####
dat$site <- factor(dat$site, levels = c("NHC", "PM", "CBP", "WB", "WBP","UNHC"))
gas <- dat  %>%
  # mutate(date = as.Date(group))%>%
  mutate(O2.ugL = DO.obs * 1000) %>%
  select(date, datetime, site, CH4.ugL, CO2.ugL, N2O.ugL, O2.ugL) %>%
  pivot_longer(ends_with("ugL"), names_to = "gas", names_pattern ='([0-9A-Z]+).',
               values_to = "concentration_ugl")

flux <- dat %>%
  filter(!is.na(datetime),
         site != "MC751") %>%
  # mutate(date = as.Date(group))%>%
  select(datetime, group,date, site,  habitat, distance_upstream_m, watertemp_C,
         depth, GPP, ER, discharge, DO.obs, no3n_mgl, ends_with('ugld')) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas",
               names_pattern ='([0-9A-Z]+).', values_to = "flux_ugld") %>%
  left_join(gas, by = c("date",'datetime', "site", "gas")) %>%
  mutate(gas = factor(gas, levels = c('CO2','O2',  'CH4', 'N2O')),
         flux_mgm2d = flux_ugld * depth)%>%
  arrange(date, distance_upstream_m)


gas_do <- dat %>%
  filter(!is.na(datetime),
         site != "MC751") %>%
  # mutate(date = as.Date(group))%>%
  mutate(DO.ugL = DO.obs * 1000) %>%
  select(datetime, date, site,  habitat, distance_upstream_m, watertemp_C,
         GPP, ER, discharge, ends_with('.ugL')) %>%
  pivot_longer(ends_with("ugL"), names_to = "gas", names_pattern ='([0-9A-Z]+).',
               values_to = "concentration_ugl") %>%
  mutate(gas = factor(gas, levels = c('DO', 'CH4', 'CO2', 'N2O')))

ggplot(flux, aes(x = -distance_upstream_m, y = no3n_mgl,
                 color = factor(date), group = factor(date))) +
  geom_line(lwd = 1.2) +
  geom_point() +
  xlab("Distance (m, upstream -> downstream)")+
  theme_bw()

png("figures/ghgconc_longitudinal_bydate_color.png", height = 5, width = 7,
    units = "in", res = 300)
    ggplot(flux, aes(x = -distance_upstream_m, y = concentration_ugl,
                     color = factor(date), group = factor(date))) +
        geom_line(lwd = 1) +
        geom_point() +
        xlab("Distance (m, upstream -> downstream)")+
        facet_wrap(gas~., scales = "free_y", ncol = 1, strip.position = "right") +
        theme_bw() +
        ylab("Concentration (ug/L)")
dev.off()

png("figures/ghgflux_longitudinal_bydate.png", height = 5, width = 7,
    units = "in", res = 300)
    ggplot(flux, aes(x = -distance_upstream_m, y = flux_ugld,
                     color = factor(date))) +
        geom_line() +
        geom_point() +
        xlab("Distance (m, upstream -> downstream)")+
        facet_wrap(gas~., scales = "free_y", ncol = 1)
dev.off()
#dumb stuff ####

# gg <- flux %>% filter(gas == 'N2O') %>%
#   select(-concentration_ugl) %>%
#   rename(concentration_ugl = flux_mgm2d)
# abs_max <- max(gg$concentration_ugl, na.rm = T)
# # get the highest point for each class
# maxs <- gg %>%
#   group_by(site) %>%
#   summarise(concentration_ugl = max(concentration_ugl) + 0.05 * abs_max)
# # get Tukey HSD results
# Tukey_test <- aov(concentration_ugl~site, data=gg) %>%
#   HSD.test("site", group=TRUE) %>%
#   .$groups %>%
#   as_tibble(rownames="site") %>%
#   rename("Letters_Tukey"="groups") %>%
#   mutate(concentration_ugl = abs_max * 1.05)
# # plot
# ggplot(gg) +
#   aes(x=site, y=concentration_ugl) +
#   geom_boxplot() +
#   geom_text(data=Tukey_test, aes(label=Letters_Tukey))

# panels for Figure 1 ####
labframe <- data.frame(gas = c('CO2', 'CH4', 'N2O', 'O2'),
         l1 = c('CO', 'CH', 'N', 'O'),
         ln = c(2, 4, 2, 2),
         l2 = c('', '', 'O', '')) %>%
    mutate(gas = factor(gas, levels = c('CO2', 'O2', 'CH4', 'N2O')))

display.brewer.all()
RColorBrewer::display.brewer.pal(8, 'Set2')
pal <- brewer.pal(8, 'Set2')

tiff("figures/final/ghgboxplots_longitudinal_boxplots.tiff", height = 4, width = 6.2,
    units = "in", res = 300)
ff <- flux %>%
   mutate(flux_mgm2d =
            ifelse(gas %in% c('CH4', 'N2O'), flux_mgm2d,
                   flux_mgm2d/1000)) %>%
    left_join(labframe, by = 'gas') %>%
  ggplot(aes(x = site, y = flux_mgm2d, fill = gas,
             group = interaction(gas, site)))+
    geom_hline(yintercept = 0, lwd = .02)+
    geom_boxplot( position = 'identity') +
    scale_fill_manual(values = c(alpha(pal[2]), alpha(pal[3]),
                                 alpha(pal[4]), alpha(pal[5])))+
    scale_x_discrete(limits = rev(levels(flux$site)),
                     labels = c(0, 2330, 2500, 5000, 6880, 8450))+
    labs(x = "Sites, distance downstream (m)",
         y = expression(paste('mg ', m^-2, ' ', d^-1,
                              '                            g ',
                              m^-2, ' ', d^-1)))+
    facet_grid(gas + l1 + ln + l2~., scales = "free_y",
               labeller = label_bquote(rows = .(l1)[.(ln)]~.(l2))) +
    ggtitle('Gas Flux Rates') +
    theme_bw() +
    theme( plot.margin = unit(c(0,0,0,.6), "lines"),
           plot.title = element_text(size=10),
           axis.title = element_text(size = 10),
           legend.position = 'none')
cc <- flux %>%
   mutate(concentration_ugl =
            ifelse(gas %in% c('CH4', 'N2O'), concentration_ugl,
                   concentration_ugl/1000)) %>%
    left_join(labframe, by = 'gas') %>%
  ggplot( aes(x = site, y = concentration_ugl, fill = gas,
              group = interaction(gas, site)))+
    geom_boxplot(position = 'identity') +
    scale_fill_manual(values = c(alpha(pal[2]), alpha(pal[3]),
                                 alpha(pal[4]), alpha(pal[5])))+
    scale_x_discrete(limits = rev(levels(flux$site)), labels = c(0, 2330, 2500, 5000, 6880, 8450))+
    labs(x = "Sites, distance downstream (m)",
         y = expression(paste(mu, 'g ', L^-1,
                              '                                  mg ', L^-1))) +
    facet_grid(gas + l1 + ln + l2~., scales = "free_y",
               labeller = label_bquote(rows = .(l1)[.(ln)]~.(l2))) +
    ggtitle('Gas Concentrations') +
    theme_bw() +
    theme( plot.margin = unit(c(0,0,0,.6), "lines"),
           plot.title = element_text(size=10),
           axis.title = element_text(size = 10),
           legend.position = 'none')
ggarrange(cc, ff, ncol = 2, align = 'h', labels = c('a', 'b'),
          label.x = 0.03, label.y = 1.02)
dev.off()


flux %>%
  ggplot(aes(x = group, y = no3n_mgl, group = factor(group)))+
  geom_hline(yintercept = 0, lwd = .2)+
  geom_boxplot(fill = 'gray70', position = 'identity') +
  theme_bw()

tiff('figures/final/CO2flux_equivalents.tif', width = 5, height = 5,
     units = 'in', res = 300)
    flux %>%
        mutate(flux.equivs = case_when(gas == 'CO2' ~ flux_mgm2d/1000,
                                       gas == 'N2O' ~ flux_mgm2d *298/1000,
                                       gas == 'CH4' ~ flux_mgm2d * 25/1000,
                                       TRUE ~ flux_mgm2d)) %>%
        filter(gas != 'O2') %>%
        group_by(group, gas) %>%
        summarize(CO2_flux_equivalents = median(flux.equivs)) %>%
        ggplot(aes(group, y = CO2_flux_equivalents, fill = gas)) +
        geom_bar(stat = 'identity') +
        # geom_text(aes(label = round(CO2_flux_equivalents, 1)), vjust=0) +
        scale_fill_brewer(type = 'qual', palette = 7)+
        # scale_fill_manual(values = c(alpha(pal[2], 0.7), alpha(pal[4], 0.7),
        #                              alpha(pal[5], 0.7))) +
        theme_minimal() +
        xlab('')+
        ylab(expression(paste(CO[2], ' flux equivalents (g  ', m^-2, d^-1, ')')))+
        geom_hline(yintercept = 0, size = .7)
dev.off()

    flux %>%
        mutate(flux.equivs = case_when(gas == 'CO2' ~ flux_mgm2d/1000,
                                       gas == 'N2O' ~ flux_mgm2d *298/1000,
                                       gas == 'CH4' ~ flux_mgm2d * 25/1000,
                                       TRUE ~ flux_mgm2d)) %>%
        filter(gas != 'O2') %>%
        group_by(group, gas) %>%
        summarize(CO2_flux_equivalents = median(flux.equivs)) %>%
        mutate(sign = case_when(CO2_flux_equivalents <0 ~ 'neg',
                                TRUE ~ 'pos')) %>%
        group_by(group, sign) %>%
        summarize(tot = round(sum(CO2_flux_equivalents), 2))

    0.35/4.33


tiff("figures/final/ghgconc_longitudinal_boxplots.tif", height = 5, width = 3.2,
    units = "in", res = 300)
    no <- flux %>%
      group_by(group) %>%
      summarize(no3_med = median(no3n_mgl, na.rm = T))%>%
      ggplot(aes(x = group, no3_med))+
      geom_line() +
      # summarize(no3_h = quantile(no3n_mgl, 0.75, na.rm = T),
      #           no3_l = quantile(no3n_mgl, 0.25, na.rm = T))%>%
      # ggplot(aes(x = group, no3_h))+
      # geom_ribbon(aes(ymin = no3_l, ymax = no3_h), fill = 'grey') +
      # ggplot(aes(x = group, no3n_mgl, group = factor(group)))+
      geom_point(data = flux, aes(group, no3n_mgl)) +
      labs(y = expression(paste(NO[3], '-N (mg ', L^-1, ')')),
           x = '', size = 9.5)+
      # geom_violin(fill = 'gray70', position = 'identity') +
      theme_bw()
    nn <- flux %>%
      filter(gas == 'N2O')%>%
       # mutate(flux_mgm2d =
       #          ifelse(gas %in% c('CH4', 'N2O'), flux_mgm2d,
       #                 flux_mgm2d/1000)) %>%
      ggplot( aes(x = group, y = flux_mgm2d, group = factor(group)))+
        geom_boxplot(fill = 'gray70', position = 'identity') +
        geom_hline(yintercept = 0, lwd = .5) +
        labs(x="Date",
               y = expression(paste(N[2],'O flux (mg ', m^-2, ' ', d^-1, ')')),
             size = 9.5) +
        theme_bw()
    ggarrange(no, nn, ncol = 1, heights = c(.5, 1), labels = c('a', 'b'),
              label.y = 1.05)
dev.off()
flux %>%
  mutate(flux_mgm2d =
           ifelse(gas %in% c('CH4', 'N2O'), flux_mgm2d,
                  flux_mgm2d/1000)) %>%
  ggplot( aes(x = group, y = flux_mgm2d, group = factor(group)))+
  geom_boxplot(fill = 'gray70', position = 'identity') +
  xlab("Date")+
  geom_hline(yintercept = 0, lwd = .5) +
  ylab('N2O flux (mg/m2/d)') +
  facet_wrap(gas~., strip.position = "right", ncol = 1, scales = 'free_y') +
  theme_bw()
png("figures/ghgflux_longitudinal_boxplots.png", height = 5.15, width = 3.2,
    units = "in", res = 300)
 CO <- flux %>%
   mutate(lab = 'CO2                      O2') %>%
   filter(gas %in% c('CO2', 'O2')) %>%
  ggplot( aes(x = site, y = flux_mgm2d/1000, group = interaction(gas, site)))+
    geom_boxplot(fill = 'gray70', position = 'identity') +
    scale_x_discrete(limits = rev(levels(flux$site)))+
    ylim(-10,10)+    ylab('g/m2/d')+
    geom_hline(yintercept = 0)+
    ggtitle('Gas Flux Rates') +
    facet_wrap(lab~., strip.position = "right", labeller = labeller(all = "CO2 O2")) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1,0,0,0), "lines"))
 CN <- flux %>%
   filter(gas %in% c('CH4', 'N2O')) %>%
  ggplot( aes(x = site, y = flux_mgm2d, group = interaction(gas, site)))+
    geom_boxplot(fill = 'gray70', position = 'identity') +
    scale_x_discrete(limits = rev(levels(flux$site)), labels = c(0, 2330, 2500, 5000, 6880, 8450))+
    xlab("Sites, distance downstream (m)")+
    ylab('mg/m2/d')+
    facet_wrap(gas~., scales = "free_y", ncol = 1, strip.position = "right") +
    theme( plot.margin = unit(c(0,0,0,0), "lines"))+
    theme_bw()
ggarrange(CO, CN, ncol = 1, align = 'v', label.y = 'Gas Flux Rate')
dev.off()


# Temporal gas and drivers Figure 3 ####
  # ggplot(flux, aes(x = as.factor(date), y = concentration_ugl)) +
  #   geom_boxplot() +
  #   facet_wrap(gas~., scales = "free_y", ncol = 1, strip.position = 'right') +
  #   theme_bw()
pds <- dvs %>%
    select(site, date, GPP, ER, discharge, watertemp_C) %>%
    filter(site != 'MC751') %>%
    group_by(date) %>%
    summarize(discharge = discharge[which(site == 'NHC')],
              watertemp_C = mean(watertemp_C, na.rm = T),
              across(any_of(c('GPP', 'ER')),
                     .fns = list(mean = ~mean(., na.rm = T),
                                 c2.5 = ~quantile(., 0.025, na.rm = T),
                                 c97.5 = ~quantile(., .975, na.rm = T)),
                     .names = '{col}_{fn}'))

pg <- gas_do %>%
    select(site, date, gas, concentration_ugl) %>%
    pivot_wider(values_from = 'concentration_ugl', names_from = 'gas') %>%
    mutate(DO = DO/1000,
           CO2 = CO2/1000,
           ch4co2 = CH4/CO2,
           o2co2 = DO/CO2) %>%
    left_join(select(dvs, site, date, doc_mgl, no3n_mgl)) %>%
    mutate(group = case_when(date == as.Date('2019-12-04') ~ as.Date('2019-12-03'),
                             date == as.Date('2020-01-30') ~ as.Date('2020-01-29'),
                             TRUE ~ date)) %>%
    group_by(group)%>%
    mutate(across(-c(site, date),
                  .fns = list(ci25 = ~quantile(., .25, na.rm = T),
                              ci75 = ~quantile(., .75, na.rm = T),
                              med = ~median(., na.rm = T)),
                  .names = '{col}_{fn}')) %>%
    arrange(date)


ll <- left_join(pds, pg) %>%
   mutate(doy = as.numeric(format(date, '%j'))) %>%
   summarize(across(-date, ~na.approx(., doy)))
 loess_mod <- loess(CO2 ~ doy, data = ll, span = 0.02)
 # ll$co2_l <- predict(loess_mod, newdata = ll)

library(viridis)
col.ER = plasma(7)[5]
col.GPP = 'forestgreen'

png('figures/figure3_timeseries_gas_conc_drivers_small.png', width = 3.5,
     height = 4.9, res = 300, units = 'in', family = 'cairo')
    m <- matrix(c(1,1,2,3,4,5,6,7), ncol = 1)
    layout(m)
    par(mar = c(.3,0,0,0), oma = c(4, 5, 1, 1))
    plot(pds$date, pds$GPP_mean, type = 'l', lwd = 2, col = 'grey35',
      ylim = c(-8.5,3.5), cex.axis = 0.8, xaxt = 'n')
    lines(pds$date, -pds$ER_mean, lwd = 2, col = 'grey35')
    polygon(c(pds$date, rev(pds$date)),
         na.approx(rollmean(c(pds$GPP_c2.5, rev(pds$GPP_c97.5)), 3, na.pad = T),
                   na.rm = F),
         border = NA, col = alpha(col.GPP, .4))
    polygon(c(pds$date, rev(pds$date)),
         na.approx(rollmean(c(-pds$ER_c2.5, rev(-pds$ER_c97.5)), 3, na.pad = T),
                   na.rm = F),
         border = NA, col = alpha(col.ER, .4))
    abline(h = 0)
    mtext('Metabolism', 2,3.2,cex = 0.63)
    mtext(expression(paste("(g"~O[2]~"m"^"-2"*" y"^"-1"*')')), 2, 2, cex = .63)
    legend('bottomright', legend = c('GPP', 'ER'),cex = .8,
        col = 'grey35', inset = .065, lty = 1, bty = 'n',  lwd = 2)
    polygon(pds$date[c(102,113,113,102)], c(-5,-5,-6,-6),
         col = alpha(col.GPP, .4), border = NA)
    polygon(pds$date[c(102,113,113,102)], c(-7.34,-7.34,-6.34,-6.34),
         col = alpha(col.ER, .4), border = NA)
    plot(pg$date, pg$CO2, pch = 20, xaxt = 'n', xlab = '', ylab = '',
      ylim = c(-0.2,8), cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$CO2_ci25, rev(pg$CO2_ci75)),
         col = alpha('grey', .5), border = NA)
    mtext(expression(paste('C'*O[2])), 2,3.2, cex = .63)
    mtext(expression(paste('(mg l'^'-1'*')')), 2,2, cex = .63)
    plot(pg$date, pg$DO, pch = 20, xaxt = 'n', xlab = '', ylab = '',
      ylim = c(6,13), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$DO_ci25, rev(pg$DO_ci75)),
         col = alpha('grey', .5), border = NA)
    axis(2, at = c(7,9,11), cex.axis = 0.8)
    mtext(expression(O[2]), 2,3.2, cex = .63)
    mtext(expression(paste('(mg l'^'-1'*')')), 2,2, cex = .63)
    plot(pg$date, pg$CH4, pch = 20, xaxt = 'n', xlab = '', ylab = '',
      ylim = c(-1, 22), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$CH4_ci25, rev(pg$CH4_ci75)),
         col = alpha('grey', .5), border = NA)
    axis(2, at = c(0,10,20), cex.axis = 0.8)
    mtext(expression(paste('C'*H[4])), 2,3.2, cex = .63)
    mtext(expression(paste('('*mu*'g l'^'-1'*')')), 2,2, cex = .63)
    plot(pg$date, pg$N2O, pch = 20, xaxt = 'n', xlab = '', ylab = '',
      ylim = c(-.1, .93), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$N2O_ci25, rev(pg$N2O_ci75)),
         col = alpha('grey', .5), border = NA)
    axis(2, at = c(0, 0.4, 0.8), cex.axis = 0.8)
    mtext(expression(paste(N[2]*'O')), 2,3.2, cex = .63)
    mtext(expression(paste('('*mu*'g l'^'-1'*')')), 2,2, cex = .63)
    plot(pds$date, pds$watertemp_C, type = 'l', lwd = 1.2, col = 'grey 25',
      xaxt = 'n', xlab = '', ylab = '', cex.axis = 0.8)
    mtext('Temp', 2,3.2, cex = .63)
    mtext(expression(paste(degree, 'C')), 2,2, cex = .63)
    plot(pds$date, pds$discharge, type = 'l', lwd = 1.2, col = 'grey25',
      xaxt = 'n', ylab = '', xlab = '', log = 'y', yaxt = 'n', cex.axis = 0.8)
    mtext('Discharge', 2,3, cex = .63)
    mtext(expression(paste("(m"^"3"*"/s)")), 2,2, cex = .63)
    axis(1, at = seq(as.Date('2019-12-01'), by = 'month', length.out = 4),
      labels = c('Dec-2019', 'Jan-2020', 'Feb-2020', 'Mar-2020'), cex.axis = 0.8)
    axis(2, at = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100), cex.axis = 0.8)
dev.off()

 # Version two of this figure, where the upstream and downstream temperature
 # and discharge are plotted separately
 # update this after the above one is super finzlized
tiff('figures/final/figure3_timeseries_gas_conc_drivers_small_v2.tif', width = 3.5,
    height = 4.9, res = 300, units = 'in', family = 'cairo')
    m <- matrix(c(1,1,2,3,4,5,6,7), ncol = 1)
    layout(m)
    par(mar = c(.3,0,0,0), oma = c(4, 4.5, 1, 1))
    plot(pds$date, pds$GPP_mean, type = 'l', lwd = 2, col = 'grey35',
         ylim = c(-8.5,3.5), cex.axis = 0.8, xaxt = 'n', xlab = '')
    lines(pds$date, -pds$ER_mean, lwd = 2, col = 'grey35')
    polygon(c(pds$date, rev(pds$date)),
            na.approx(rollmean(c(pds$GPP_c2.5, rev(pds$GPP_c97.5)), 3, na.pad = T),
                      na.rm = F),
            border = NA, col = alpha(col.GPP, .4))
    polygon(c(pds$date, rev(pds$date)),
            na.approx(rollmean(c(-pds$ER_c2.5, rev(-pds$ER_c97.5)), 3, na.pad = T),
                      na.rm = F),
            border = NA, col = alpha(col.ER, .4))
    abline(h = 0)
    mtext('Metabolism', 2,3.2,cex = 0.63)
    mtext(expression(paste("(g"~O[2]~m^-2~y^-1,")")), 2, 2, cex = .63)
    legend('bottomright', legend = c('GPP', 'ER'),cex = .8,
           col = 'grey35', inset = .065, lty = 1, bty = 'n',  lwd = 2)
    polygon(pds$date[c(103,114,114,103)], c(-5,-5,-6,-6),
            col = alpha(col.GPP, .4), border = NA)
    polygon(pds$date[c(103,114,114,103)], c(-7.34,-7.34,-6.34,-6.34),
            col = alpha(col.ER, .4), border = NA)
    plot(pg$date, pg$CO2, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-0.2,8), cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$CO2_ci25, rev(pg$CO2_ci75)),
            col = alpha('grey', .5), border = NA)
    mtext(expression(CO[2]), 2, 3.2, cex = .63)
    mtext(expression(paste('(mg'~L^-1,')')), 2, 2, cex = .63)
    plot(pg$date, pg$DO, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(6,13), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$DO_ci25, rev(pg$DO_ci75)),
            col = alpha('grey', .5), border = NA)
    axis(2, at = c(7,9,11), cex.axis = 0.8)
    mtext(expression(O[2]), 2, 3.2, cex = .63)
    mtext(expression(paste('(mg'~L^-1,')')), 2, 2, cex = .63)
    plot(pg$date, pg$CH4, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-1, 22), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$CH4_ci25, rev(pg$CH4_ci75)),
            col = alpha('grey', .5), border = NA)
    axis(2, at = c(0,10,20), cex.axis = 0.8)
    mtext(expression(CH[4]), 2, 3.2, cex = .63)
    mtext(expression(paste('('*mu*'g'~L^-1,')')), 2,2, cex = .63)
    plot(pg$date, pg$N2O, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-.1, .93), yaxt = 'n', cex.axis = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$N2O_ci25, rev(pg$N2O_ci75)),
            col = alpha('grey', .5), border = NA)
    axis(2, at = c(0, 0.4, 0.8), cex.axis = 0.8)
    mtext(expression(paste(N[2]*'O')), 2, 3.2, cex = .63)
    mtext(expression(paste('('*mu*'g'~L^-1,')')), 2,2, cex = .63)
    plot(dvs[dvs$site == "UNHC",]$date, dvs[dvs$site == 'UNHC',]$watertemp_C,
         type = 'l', lwd = 1.2, col = 'grey 10',
         xaxt = 'n', xlab = '', ylab = '', cex.axis = 0.8)
    lines(dvs[dvs$site == "NHC",]$date, dvs[dvs$site == 'NHC',]$watertemp_C,
         lwd = 1.2, col = 'grey 10', lty = 2)
    mtext("Temp", 2, 3.2, cex = .63)
    mtext(expression(paste('('~degree, 'C)')), 2,2, cex = .63)
    legend('topleft', c('upstream', 'downstream'), lty = c(1,2),
           col = 'grey 10', bty = 'n', cex = .7)
    plot(dvs[dvs$site == "UNHC",]$date, dvs[dvs$site == 'UNHC',]$discharge,
         type = 'l', lwd = 1.2, col = 'grey 10', log = 'y', xaxt = 'n',
         yaxt = 'n', xlab = '', ylab = '', cex.axis = 0.8, ylim = c(0.1, 500))
    legend('topleft', c('upstream', 'downstream'), lty = c(1,2),
           col = 'grey 10', bty = 'n', cex = .7)
    lines(dvs[dvs$site == "NHC",]$date, dvs[dvs$site == 'NHC', ]$discharge,
         lwd = 1.2, col = 'grey 10', lty = 2)
    mtext("Discharge", 2, 3.2, cex = .63)
    mtext(expression(paste("(m"^"3"~s^-1,")")), 2,2, cex = .63)
    axis(1, at = seq(as.Date('2019-12-01'), by = 'month', length.out = 4),
         labels = c('Dec-2019', 'Jan-2020', 'Feb-2020', 'Mar-2020'), cex.axis = 0.8)
    axis(2, at = c(0.1, 1, 100), labels = c(0.1, 1, 100), cex.axis = 0.8)
dev.off()

 # Version three of this figure, with some panels combined per reviewer request
ratios <- read_csv('data/fraction_of_instream_production_CO2_and_CH4CO2ratios.csv') %>%
    mutate(DO_extra = DO.obs - DO.sat,
           O2CO2_extra = DO_extra/CO2.umol_extra)%>%
    filter(CO2.umol_extra != 0) %>%
    group_by(group)%>%
    summarize(ch4co2 = mean(CH4CO2_extra, na.rm = T),
              o2co2 = mean(O2CO2_extra, na.rm = T))



tiff('figures/final/figure3_timeseries_gas_conc_drivers_small_v3.tif', width = 3.5,
    height = 5.2, res = 300, units = 'in', family = 'cairo')
    m <- matrix(c(1,1,2,2,3,3,4,4,5,5), ncol = 1)
    layout(m)
    par(mar = c(.3,0,0,0), oma = c(5, 4, 1, 3))
    plot(pds$date, pds$GPP_mean, type = 'l', lwd = 2, col = 'grey35',
         ylim = c(-8.5,3.5), cex.axis = 0.8, xaxt = 'n', xlab = '', padj = 0.8)
    lines(pds$date, -pds$ER_mean, lwd = 2, col = 'grey35')
    polygon(c(pds$date, rev(pds$date)),
            na.approx(rollmean(c(pds$GPP_c2.5, rev(pds$GPP_c97.5)), 3, na.pad = T),
                      na.rm = F),
            border = NA, col = alpha(pal[1], .7))
    polygon(c(pds$date, rev(pds$date)),
            na.approx(rollmean(c(-pds$ER_c2.5, rev(-pds$ER_c97.5)), 3, na.pad = T),
                      na.rm = F),
            border = NA, col = alpha(pal[7], .7))
    abline(h = 0)
    mtext('Metabolism', 2,2.7,cex = 0.63)
    mtext(expression(paste("(g"~O[2]~m^-2~y^-1,")")), 2, 1.5, cex = .63)
    legend('bottomright', legend = c('GPP', 'ER'),cex = .8,
           col = 'grey35', inset = .065, lty = 1, bty = 'n',  lwd = 2)
    polygon(pds$date[c(101,112,112,101)], c(-4.5,-4.5,-5.5,-5.5),
            col = alpha(pal[1], .7), border = NA)
    polygon(pds$date[c(101,112,112,101)], c(-7,-7,-6,-6),
            col = alpha(pal[7], .7), border = NA)

    plot(dvs[dvs$site == "UNHC",]$date, dvs[dvs$site == 'UNHC',]$watertemp_C,
         type = 'l', lwd = 1.5, col = 'grey40', ylim = c(-18,24),
         xaxt = 'n', xlab = '', ylab = '', cex.axis = 0.8, yaxt = 'n')
    lines(dvs[dvs$site == "NHC",]$date, dvs[dvs$site == 'NHC',]$watertemp_C,
         lwd = 1.5, col = 'grey40', lty = 2)
    mtext(expression(paste('Temp ('~degree, 'C)')), 4, 1.6, cex = .63, adj = 1)
    axis(4, at = c(5, 10, 15, 20), cex.axis = 0.8, padj = -1)
    legend('topleft', c('upstream', 'downstream'), lty = c(1,2),
           col = 'grey40', bty = 'n', cex = .7)
    par(new = T)
    plot(dvs[dvs$site == "UNHC",]$date, dvs[dvs$site == 'UNHC',]$discharge,
         type = 'l', lwd = 1.5, col = 'grey40', log = 'y', xaxt = 'n',
         yaxt = 'n', xlab = '', ylab = '', cex.axis = 0.8, ylim = c(.2, 150))
    lines(dvs[dvs$site == "NHC",]$date, dvs[dvs$site == 'NHC', ]$discharge,
         lwd = 1.5, col = 'grey40', lty = 2)
    mtext("Discharge", 2, 2.7, cex = .63, adj = 0.01)
    mtext(expression(paste("(m"^"3"~s^-1,")")), 2, 1.5, cex = .63, adj = 0.1)
    axis(2, at = c(0.2, 1,5), labels = c(0.2, 1, 5), cex.axis = 0.8, padj = 0.8)
    plot(pg$date, pg$CO2, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-0.2,8.9), cex.axis = 0.8, col = pal[2], padj = 0.8)
    polygon(c(pg$date, rev(pg$date)), c(pg$CO2_ci25, rev(pg$CO2_ci75)),
            col = alpha(pal[2], .5), border = NA)
    mtext(expression(paste(CO[2], ' (mg'~L^-1,')')), 2, 1.5, cex = .63)
    par(new = T)
    plot(pg$date, pg$DO, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(6,12.5), yaxt = 'n', cex.axis = 0.8, col = pal[3])
    polygon(c(pg$date, rev(pg$date)), c(pg$DO_ci25, rev(pg$DO_ci75)),
            col = alpha(pal[3], .5), border = NA)
    axis(4, at = c(7,9,11), cex.axis = 0.8, padj = -1)
    mtext(expression(paste(O[2], ' (mg'~L^-1,')')), 4, 1.6, cex = .63)
    par(new = T)
    plot(ratios$group, ratios$o2co2, type = 'l', lty = 2, col = pal[3], lwd = 1.5,
         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    plot(pg$date, pg$CH4, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-1, 22), yaxt = 'n', cex.axis = 0.8, col = pal[4])
    polygon(c(pg$date, rev(pg$date)), c(pg$CH4_ci25, rev(pg$CH4_ci75)),
            col = alpha(pal[4], .5), border = NA)
    axis(2, at = c(0,10,20), cex.axis = 0.8, padj = 0.8)
    mtext(expression(paste(CH[4], ' ('*mu*'g'~L^-1,')')), 2,1.5, cex = .63)
    par(new = T)
    plot(ratios$group, ratios$ch4co2, type = 'l', lty = 2, col = pal[4], lwd = 1.5,
         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    par(new = T)
    plot(pg$date, pg$N2O, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         ylim = c(-.1, .93), yaxt = 'n', cex.axis = 0.8, col = pal[5])
    polygon(c(pg$date, rev(pg$date)), c(pg$N2O_ci25, rev(pg$N2O_ci75)),
            col = alpha(pal[5], .5), border = NA)
    axis(4, at = c(0, 0.4, 0.8), cex.axis = 0.8, padj = -1)
    mtext(expression(paste(N[2]*'O', ' ('*mu*'g'~L^-1,')')), 4,1.6, cex = .63)

    plot(pg$date, pg$doc_mgl, pch = 20, xaxt = 'n', xlab = '', ylab = '',
          yaxt = 'n', cex.axis = 0.8, col = pal[1])
    polygon(c(pg$date, rev(pg$date)), c(pg$doc_mgl_ci25, rev(pg$doc_mgl_ci75)),
            col = alpha(pal[1], .5), border = NA)
    axis(2, at = c(2,4,6, 8), cex.axis = 0.8, padj = 0.8)
    mtext(expression(paste('DOC (mg ', L^-1, ')')), 2,1.5, cex = .63)
    par(new = T)
    plot(pg$date, pg$no3n_mgl, pch = 20, xaxt = 'n', xlab = '', ylab = '',
         yaxt = 'n', cex.axis = 0.8, col = pal[8])
    polygon(c(pg$date, rev(pg$date)), c(pg$no3n_mgl_ci25, rev(pg$no3n_mgl_ci75)),
            col = alpha(pal[8], .5), border = NA)
    axis(4, at = c(0, 0.4, 0.8), cex.axis = 0.8, padj = -1)
    mtext(expression(paste(NO[3],'-N (mg ', L^-1,')')), 4,1.6, cex = .63)
    axis(1, at = seq(as.Date('2019-12-01'), by = 'month', length.out = 4),
         labels = c('Dec-2019', 'Jan-2020', 'Feb-2020', 'Mar-2020'),
         cex.axis = 0.9,  padj = -0.6)
    legend('bottomleft', inset = c(0, -0.8), xpd = NA, ncol = 4, bty = 'n',
           legend = c(expression(CO[2]), expression(O[2]),
                      expression(CH[4]), expression(N[2]*O),
                      'DOC', expression(NO[3]*-N),
                      expression(paste('excess ', O[2], ':', CO[2])),
                      expression(paste('excess ', CH[4], ':', CO[2]))),
           col = pal[c(2,3,4,5,1,8,3,4)], pch = c(rep(19,6), NA, NA),
           lty = c(rep(NA, 6), 2,2), lwd = 1.5, x.intersp = .1,cex = 1,
           text.width = c(rep(13,4), rep(14,2), 18, 18))
dev.off()


 # png("figures/ghgconc_temporal.png", height = 3.5, width = 6, units = "in",
 #    res = 300)
  ggplot(gas_do, aes(x = date, y = concentration_ugl)) +
    geom_point(show.legend = F) +
    geom_smooth(show.legend = F, col = 'black') +
    facet_wrap(gas~., scales = "free_y", ncol = 1, strip.position = 'right') +
    theme_bw()
# dev.off()
# png("figures/met_temporal.png", height = 1.8, width = 6, units = "in",
#     res = 300)
  ggplot(dvs, aes(x = date, y = GPP, group = site)) +
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = -ER), col = 'sienna') +
    geom_hline(yintercept = 0, lwd = .2) +
    ylab('Met gC/m2/d') +
    theme_bw()
  # dev.off()
# png("figures/QT_temporal.png", height = 2.2, width = 6, units = "in",
#     res = 300)
  dvs %>%
    filter(site == 'NHC') %>%
    mutate(logQ = log(discharge)) %>%
    pivot_longer(cols = c(logQ, watertemp_C), names_to = "var",
                 values_to = "val") %>%

  ggplot(aes(x = date, y = val)) +
    geom_line() +
    ylab('Water Temperature (C)                                   log(discharge)')+
    facet_wrap(~var, ncol = 1, scales = 'free_y', strip.position = 'right')+
    theme_bw()
    # dev.off()

  # mutate(date = as.Date(group))%>%
# dev.off()


ggplot(flux, aes(date, GPP)) +
  geom_line() +
  geom_point()+
  facet_wrap(~site)
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
