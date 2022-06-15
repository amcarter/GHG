# plot fraction of CO2 produced instream
setwd('C:/Users/alice.carter/git/ghg_patterns_nhc/')
# devtools::install_github("psyteachr/introdataviz")
library(tidyverse)
library(ggpubr)
library(introdataviz)

dat <- read_csv('data/fraction_of_instream_production_CO2_and_CH4CO2ratios.csv') %>%
    mutate(distance_m = 8450 - distance_upstream_m,
           distance_m = factor(distance_m))
dvs <- read_csv("data/ghg_filled_drivers_dataframe.csv") %>%
    filter(site == 'UNHC')
           date >= min(unique(dat$group)),
           date <= max(unique(dat$group)))

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
png('figures/CO2_flux_from_NEP_site_violins.png', width = 6.5, height = 4,
    res = 300,  units = 'in')
ggplot(dd, aes(distance_m, gm2d, fill = category) )+
    introdataviz::geom_split_violin(alpha = .4, trim = FALSE, adjust = 1.1) +
    introdataviz::geom_split_violin(aes(y = gm2d_low),
                                    alpha = .4, trim = FALSE, adjust = 1.1) +
    introdataviz::geom_split_violin(aes(y = gm2d_high),
                                    alpha = .4, trim = FALSE, adjust = 1.1) +
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
                                 expression(paste(CO[2], ' from NEP with RQ = 0.6, 0.8, 1                 Distribution median')),
                                 'median'))+
    labs(x = 'Distance Downstream (m)',
         y = expression(paste(CO[2], ' flux (g/', m^2, '/d)')))+
    theme_bw()+
    theme(legend.position = 'top',
          # legend.direction = 'vertical',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()
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

png('figures/CO2_flux_from_NEP_date_bar.png', width = 6.5, height = 4,
    res = 300,  units = 'in')

dp %>%
    mutate(NEP_pos = NEP_pos + NEP_neg,
        Extra = CO2_flux - NEP_pos)%>%
    pivot_longer(cols = any_of(c('NEP_pos', 'Extra')),
                 names_to = 'category', values_to = 'gm2d') %>%
    ggplot(aes(group, gm2d, fill = category) )+
    geom_bar(stat = 'identity') + #ylim(-.1,1) +
    geom_hline(yintercept = 0, size = .3)+
    ylab('CO2 flux (g/m2/d)') +
    xlab('')+
    guides(fill = F)+
    scale_fill_manual(values = c('grey70', plasma(7)[5]))+
    geom_line(data = dvs,
              aes(x = date, y = (log(discharge)+1.75)/6, fill = NULL),
              size = .5)+
    theme_bw()+
    theme(plot.margin = unit(c(.1,.2,.2,.2), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()
