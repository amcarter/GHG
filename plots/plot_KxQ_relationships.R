# calculate and plot K600xQ curves for each of the sites
# use the expected K600 value for each time point, not the measured one

library(tidyverse)

setwd("C:/Users/alice.carter/git/ghg_patterns_nhc/")
dat <- readRDS('data/site_data/met_preds_stream_metabolizer_O2.rds')

k <- dat$preds %>%
    filter(era == 'now') %>%
    dplyr::select(date, site, year, K600, discharge)
sites <- read_csv('data/site_data/NHCsite_metadata.csv') %>%
    slice(c(1:5,7)) %>%
    mutate(dist_downstream = 8450 - distance_m,
           dist_downstream = paste(as.character(dist_downstream),
                                   'meters', sep = ' ')) %>%
    dplyr::select(site = sitecode, dist_downstream)

k <- left_join(k, sites)
# plot KQ relationships for each site
png('figures/final/KxQ_curves.png', width = 6, height = 3,
    units = 'in', res = 300)
    ggplot(k, aes(log(discharge), K600)) +
        geom_point(size = .5)+
        geom_smooth( size = .7,  method = loess, se = FALSE) +
        facet_wrap(.~dist_downstream) +
        theme_bw()
dev.off()

for(i in 1:6){
    s = sites$site[i]
    dd <- filter(k, site == s)
    l <- loess(K600 ~ log(discharge), data = dd)
    predict(l, 8.5', 'NHC_' 'NHC_0'))
}
