#Calculate the fraction of CO2 fluxes that can be explained by DOC uptake vs leaf litter

library(tidyverse)
library(viridis)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")

dat <- read_csv('data/ghg_flux_complete_drivers_dataframe.csv')

ggplot(dat, aes(date, doc_mgl , col = log10(discharge)))+
    geom_point(size = 2)+
    scale_color_gradientn(colors = plasma(8)[1:7],
                          name = expression(Discharge (m^3~s^-1)),
                          na.value = "grey80",
                          breaks = c(-0.69897,-.1549, .39794),
                          labels = c(0.2,0.7,2.5)) +
    theme_bw()+
    xlab('Date')+ ylab(expression(paste('DOC (mg ', l^-1, ')')))


v_f_mmmin = 0.26    # lechate uptake velocity from Mineau et al 2016
                    # median reach-scale ambient DOC vf
v_f_mday <- v_f_mmmin/1000*60*24
group.labs <- c('Nov-11', 'Nov-20', 'Nov-26', 'Dec-03', 'Dec-12',
                      'Jan-05', 'Jan-29', 'Feb-12', 'Feb-27', 'Mar-11', 'Mar-20')
names(group.labs) = unique(dat$group)

tiff('figures/final/DOC_uptake_vs_ER.tif', width= 7, height = 4,
     res = 300, units = 'in')
dat %>%
    select(date, group, site, doc_mgl, depth, width = width_march_m,
           discharge, ER, CO2.flux_ugld, DO.obs) %>%
    mutate(DOC_gd = doc_mgl * discharge,
           DOC_uptake = doc_mgl * v_f_mday,
           CO2_flux= CO2.flux_ugld/1000 * depth * 12/44,
           Respiration = ER * 12/32,
           res_time_d = width*depth*1000/discharge/60/60/24,
           DOC_av = DOC_gd * res_time_d/width/1000,
           lab = substr(group, 6,10)) %>%

    pivot_longer(cols = c('CO2_flux', 'DOC_uptake'), names_to = 'rate',
                     values_to = 'gCm2d')%>%
    ggplot(aes(rate, gCm2d, fill = rate)) +
    geom_boxplot(position = 'identity')+
    facet_wrap(.~group, nrow = 1, strip.position = 'bottom',
               labeller = labeller(group = group.labs))+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0, 'lines'),
          legend.title = element_blank(),
          legend.position = 'right') +
    ylab(expression(paste('flux (g C ', m^-2, d^-1, ')')))
dev.off()



