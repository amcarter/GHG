# Add whitebox slope to sites dataframe

library(tidyverse)

setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")

sites <- read_csv('data/site_data/NHCsite_metadata.csv')

# load slopes from whitebox
slope <- read_csv('data/sites_whitebox_slopes.csv') %>%
  select(sitename, slope_wbx = slope) %>%
  right_join(sites, by = 'sitename') %>%
  select(colnames(sites), slope_wbx) %>%
  rename(slope_nhd = slope)

# The slope for WBP is much higher than it should be having been to the site
# I think that the raster squares must be large enough to capture the riffle
# that is only ~ 50m downstream. Because of this, I am going to use the slope
# calculated for the power cut site, located 180 m upstream and well within
# the same pool.
slope <- slope %>%
  mutate(slope_wbx = case_when(sitecode == 'WBP' ~ 
                                 slope$slope_wbx[slope$sitecode == 'PWC'],
                               TRUE ~ slope_wbx))
                               
write_csv(slope, 'data/site_data/NHCsite_metadata.csv')